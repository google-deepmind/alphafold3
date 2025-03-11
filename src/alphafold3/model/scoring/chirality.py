# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Chirality detection and comparison."""

from collections.abc import Mapping

from absl import logging
from alphafold3 import structure
from alphafold3.constants import chemical_components
from alphafold3.data.tools import rdkit_utils
import rdkit.Chem as rd_chem

_CHIRAL_ELEMENTS = frozenset({'C', 'S'})


def _find_chiral_centres(mol: rd_chem.Mol) -> dict[str, str]:
  """Find chiral centres and detect their chirality.

  Only elements listed in _CHIRAL_ELEMENTS are considered as centres.

  Args:
    mol: The molecule for which to detect chirality.

  Returns:
    Map from chiral centre atom names to identified chirality.
  """
  chiral_centres = rd_chem.FindMolChiralCenters(
      mol, force=True, includeUnassigned=False, useLegacyImplementation=True
  )
  atom_name_by_idx = {
      atom.GetIdx(): atom.GetProp('atom_name') for atom in mol.GetAtoms()
  }
  atom_chirality_by_name = {atom_name_by_idx[k]: v for k, v in chiral_centres}
  return {
      k: v
      for k, v in atom_chirality_by_name.items()
      if any(k[: len(el)].upper() == el for el in _CHIRAL_ELEMENTS)
  }


def _chiral_match(mol1: rd_chem.Mol, mol2: rd_chem.Mol) -> bool:
  """Compares chirality of two Mols. Mol1 can match a subset of mol2."""

  mol1_atom_names = {a.GetProp('atom_name') for a in mol1.GetAtoms()}
  mol2_atom_names = {a.GetProp('atom_name') for a in mol2.GetAtoms()}
  if mol1_atom_names != mol2_atom_names:
    if not mol1_atom_names.issubset(mol2_atom_names):
      raise ValueError('Mol1 atoms are not a subset of mol2 atoms.')

  mol1_chiral_centres = _find_chiral_centres(mol1)
  mol2_chiral_centres = _find_chiral_centres(mol2)
  if set(mol1_chiral_centres) != set(mol2_chiral_centres):
    if not set(mol1_chiral_centres).issubset(mol2_chiral_centres):
      return False
  chirality_matches = {
      centre_atom: chirality1 == mol2_chiral_centres[centre_atom]
      for centre_atom, chirality1 in mol1_chiral_centres.items()
      if '?' != mol2_chiral_centres[centre_atom]
  }
  return all(chirality_matches.values())


def _mol_from_ligand_struc(
    ligand_struc: structure.Structure,
    ref_mol: rd_chem.Mol,
) -> rd_chem.Mol | None:
  """Creates a Mol object from a ligand structure and reference mol."""

  if ligand_struc.num_residues(count_unresolved=True) > 1:
    raise ValueError('ligand_struc %s has more than one residue.')
  coords_by_atom_name = dict(zip(ligand_struc.atom_name, ligand_struc.coords))

  ref_mol = rdkit_utils.sanitize_mol(
      ref_mol,
      sort_alphabetically=False,
      remove_hydrogens=True,
  )

  mol = rd_chem.Mol(ref_mol)
  mol.RemoveAllConformers()

  atom_indices_to_remove = [
      a.GetIdx()
      for a in mol.GetAtoms()
      if a.GetProp('atom_name') not in coords_by_atom_name
  ]
  editable_mol = rd_chem.EditableMol(mol)
  # Remove indices from the largest to smallest, to avoid invalidating.
  for atom_idx in atom_indices_to_remove[::-1]:
    editable_mol.RemoveAtom(atom_idx)
  mol = editable_mol.GetMol()

  conformer = rd_chem.Conformer(mol.GetNumAtoms())
  for atom_idx, atom in enumerate(mol.GetAtoms()):
    atom_name = atom.GetProp('atom_name')
    coords = coords_by_atom_name[atom_name]
    conformer.SetAtomPosition(atom_idx, coords.tolist())
  mol.AddConformer(conformer)
  try:
    rd_chem.AssignStereochemistryFrom3D(mol)
  except RuntimeError as e:
    # Catch only this specific rdkit error.
    if 'Cannot normalize a zero length vector' in str(e):
      return None
    else:
      raise
  return mol


def _maybe_mol_from_ccd(res_name: str) -> rd_chem.Mol | None:
  """Creates a Mol object from CCD information if res_name is in the CCD."""
  ccd = chemical_components.cached_ccd()
  ccd_cif = ccd.get(res_name)
  if not ccd_cif:
    logging.warning('No ccd information for residue %s.', res_name)
    return None
  try:
    mol = rdkit_utils.mol_from_ccd_cif(ccd_cif, force_parse=False)
  except rdkit_utils.MolFromMmcifError as e:
    logging.warning('Failed to create mol from ccd for %s: %s', res_name, e)
    return None
  if mol is None:
    raise ValueError('Failed to create mol from ccd for %s.' % res_name)
  mol = rdkit_utils.sanitize_mol(
      mol,
      sort_alphabetically=False,
      remove_hydrogens=True,
  )
  return mol


def compare_chirality(
    test_struc: structure.Structure,
    ref_mol_by_chain: Mapping[str, rd_chem.Mol] | None = None,
) -> dict[str, bool]:
  """Compares chirality of ligands in a structure with reference molecules.

  We do not enforce that ligand atoms exactly match, only that the ligand atoms
  and chiral centres are a subset of those in ref mol.

  Args:
    test_struc: The structure for whose ligands to match chirality.
    ref_mol_by_chain: Optional dictionary mapping chain IDs to mol objects with
      conformers to compare against. If this is not provided, the comparison is
      to the corresponding ligands in the CCD if the ligand residue name is in
      the CCD.

  Returns:
    Dictionary mapping chain id to whether chirality mismatches the ref mol.
    Only single residue ligands where reference molecules are available are
    compared.
  """
  ref_mol_by_chain = ref_mol_by_chain or {}
  test_struc = test_struc.filter_to_entity_type(ligand=True)
  name = test_struc.name
  chiral_match_by_chain_id = {}
  for chain_id in test_struc.chains:
    chain_struc = test_struc.filter(chain_id=chain_id)
    # Only compare single-residue ligands.
    if chain_struc.num_residues(count_unresolved=True) > 1:
      logging.warning('%s: Chain %s has >1 residues. Skipping.', name, chain_id)
      continue
    if chain_id not in ref_mol_by_chain:
      ref_mol = _maybe_mol_from_ccd(chain_struc.res_name[0])
    else:
      ref_mol = ref_mol_by_chain[chain_id]
    if ref_mol is None:
      logging.warning(
          '%s: Ref mol is None for chain %s. Skipping.', name, chain_id
      )
      continue
    mol = _mol_from_ligand_struc(
        ligand_struc=chain_struc,
        ref_mol=ref_mol,
    )
    if mol is None:
      logging.warning(
          '%s: Failed to create mol for chain %s. Skipping.', name, chain_id
      )
      continue
    chiral_match_by_chain_id[chain_id] = _chiral_match(mol, ref_mol)
  return chiral_match_by_chain_id
