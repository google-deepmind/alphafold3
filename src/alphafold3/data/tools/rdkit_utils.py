# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Tools for calculating features for ligands."""

import collections
from collections.abc import Mapping, Sequence

from absl import logging
from alphafold3.cpp import cif_dict
import numpy as np
import rdkit.Chem as rd_chem
from rdkit.Chem import AllChem as rd_all_chem


_RDKIT_MMCIF_TO_BOND_TYPE: Mapping[str, rd_chem.BondType] = {
    'SING': rd_chem.BondType.SINGLE,
    'DOUB': rd_chem.BondType.DOUBLE,
    'TRIP': rd_chem.BondType.TRIPLE,
}

_RDKIT_BOND_TYPE_TO_MMCIF: Mapping[rd_chem.BondType, str] = {
    v: k for k, v in _RDKIT_MMCIF_TO_BOND_TYPE.items()
}

_RDKIT_BOND_STEREO_TO_MMCIF: Mapping[rd_chem.BondStereo, str] = {
    rd_chem.BondStereo.STEREONONE: 'N',
    rd_chem.BondStereo.STEREOE: 'E',
    rd_chem.BondStereo.STEREOZ: 'Z',
    rd_chem.BondStereo.STEREOCIS: 'Z',
    rd_chem.BondStereo.STEREOTRANS: 'E',
}


class MolFromMmcifError(Exception):
  """Raised when conversion from mmCIF to RDKit Mol fails."""


class UnsupportedMolBondError(Exception):
  """Raised when we try to handle unsupported RDKit bonds."""


def _populate_atoms_in_mol(
    mol: rd_chem.Mol,
    atom_names: Sequence[str],
    atom_types: Sequence[str],
    atom_charges: Sequence[int],
    implicit_hydrogens: bool,
    ligand_name: str,
    atom_leaving_flags: Sequence[str],
):
  """Populate the atoms of a Mol given atom features.

  Args:
    mol: Mol object.
    atom_names: Names of the atoms.
    atom_types: Types of the atoms.
    atom_charges: Charges of the atoms.
    implicit_hydrogens: Whether to mark the atoms to allow implicit Hs.
    ligand_name: Name of the ligand which the atoms are in.
    atom_leaving_flags: Whether the atom is possibly a leaving atom. Values from
      the CCD column `_chem_comp_atom.pdbx_leaving_atom_flag`. The expected
      values are 'Y' (yes), 'N' (no), '?' (unknown/unset, interpreted as no).

  Raises:
    ValueError: If atom type is invalid.
  """
  # Map atom names to the position they will take in the rdkit molecule.
  atom_name_to_idx = {name: i for i, name in enumerate(atom_names)}

  for atom_name, atom_type, atom_charge, atom_leaving_flag in zip(
      atom_names, atom_types, atom_charges, atom_leaving_flags, strict=True
  ):
    try:
      if atom_type == 'X':
        atom_type = '*'
      atom = rd_chem.Atom(atom_type)
    except RuntimeError as e:
      raise ValueError(f'Failed to use atom type: {str(e)}') from e

    if not implicit_hydrogens:
      atom.SetNoImplicit(True)

    atom.SetProp('atom_name', atom_name)
    atom.SetProp('atom_leaving_flag', atom_leaving_flag)
    atom.SetFormalCharge(atom_charge)
    residue_info = rd_chem.AtomPDBResidueInfo()
    residue_info.SetName(_format_atom_name(atom_name, atom_type))
    residue_info.SetIsHeteroAtom(True)
    residue_info.SetResidueName(ligand_name)
    residue_info.SetResidueNumber(1)
    atom.SetPDBResidueInfo(residue_info)
    atom_index = mol.AddAtom(atom)
    assert atom_index == atom_name_to_idx[atom_name]


def _populate_bonds_in_mol(
    mol: rd_chem.Mol,
    atom_names: Sequence[str],
    bond_begins: Sequence[str],
    bond_ends: Sequence[str],
    bond_orders: Sequence[str],
    bond_is_aromatics: Sequence[bool],
):
  """Populate the bonds of a Mol given bond features.

  Args:
    mol: Mol object.
    atom_names: Names of atoms in the molecule.
    bond_begins: Names of atoms at the beginning of the bond.
    bond_ends: Names of atoms at the end of the bond.
    bond_orders: What order the bonds are.
    bond_is_aromatics: Whether the bonds are aromatic.
  """
  atom_name_to_idx = {name: i for i, name in enumerate(atom_names)}
  for begin, end, bond_type, is_aromatic in zip(
      bond_begins, bond_ends, bond_orders, bond_is_aromatics, strict=True
  ):
    begin_name, end_name = atom_name_to_idx[begin], atom_name_to_idx[end]
    bond_idx = mol.AddBond(begin_name, end_name, bond_type)
    mol.GetBondWithIdx(bond_idx - 1).SetIsAromatic(is_aromatic)


def sanitize_mol(mol, sort_alphabetically, remove_hydrogens) -> rd_chem.Mol:
  # https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.SanitizeMol
  # Kekulize, check valencies, set aromaticity, conjugation and hybridization.
  # This can repair e.g. incorrect aromatic flags.
  rd_chem.SanitizeMol(mol)
  if sort_alphabetically:
    mol = sort_atoms_by_name(mol)
  if remove_hydrogens:
    mol = rd_chem.RemoveHs(mol)
  return mol


def _add_conformer_to_mol(mol, conformer, force_parse) -> rd_chem.Mol:
  # Create conformer and use it to assign stereochemistry.
  if conformer is not None:
    try:
      mol.AddConformer(conformer)
      rd_chem.AssignStereochemistryFrom3D(mol)
    except ValueError as e:
      logging.warning('Failed to parse conformer: %s', e)
      if not force_parse:
        raise


def mol_from_ccd_cif(
    mol_cif: cif_dict.CifDict,
    *,
    force_parse: bool = False,
    sort_alphabetically: bool = True,
    remove_hydrogens: bool = True,
    implicit_hydrogens: bool = False,
) -> rd_chem.Mol:
  """Creates an rdkit Mol object from a CCD mmcif data block.

  The atoms are renumbered so that their names are in alphabetical order and
  these names are placed on the atoms under property 'atom_name'.
  Only hydrogens which are not required to define the molecule are removed.
  For example, hydrogens that define stereochemistry around a double bond are
  retained.
  See this link for more details.
  https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RemoveHs

  Args:
     mol_cif: An mmcif object representing a molecule.
     force_parse: If True, assumes missing aromatic flags are false, substitutes
       deuterium for hydrogen, assumes missing charges are 0 and ignores missing
       conformer / stereochemistry information.
     sort_alphabetically: True: sort atom alphabetically; False: keep CCD order
     remove_hydrogens: if True, remove non-important hydrogens
     implicit_hydrogens: Sets a marker on the atom that allows implicit Hs.

  Returns:
     An rdkit molecule, with the atoms sorted by name.

  Raises:
    MolToMmcifError: If conversion from mmcif to rdkit Mol fails. More detailed
      error is available as this error's cause.
  """
  # Read data fields.
  try:
    atom_names, atom_types, atom_charges, atom_leaving_flags = parse_atom_data(
        mol_cif, force_parse
    )
    bond_begins, bond_ends, bond_orders, bond_is_aromatics = parse_bond_data(
        mol_cif, force_parse
    )
    lig_name = mol_cif['_chem_comp.id'][0].rjust(3)
  except (KeyError, ValueError) as e:
    raise MolFromMmcifError from e

  # Build Rdkit molecule.
  mol = rd_chem.RWMol()

  # Per atom features.
  try:
    _populate_atoms_in_mol(
        mol=mol,
        atom_names=atom_names,
        atom_types=atom_types,
        atom_charges=atom_charges,
        implicit_hydrogens=implicit_hydrogens,
        ligand_name=lig_name,
        atom_leaving_flags=atom_leaving_flags,
    )
  except (ValueError, RuntimeError) as e:
    raise MolFromMmcifError from e

  _populate_bonds_in_mol(
      mol, atom_names, bond_begins, bond_ends, bond_orders, bond_is_aromatics
  )

  try:
    conformer = _parse_ideal_conformer(mol_cif)
  except (KeyError, ValueError) as e:
    logging.warning('Failed to parse ideal conformer: %s', e)
    if not force_parse:
      raise MolFromMmcifError from e
    conformer = None

  mol.UpdatePropertyCache(strict=False)

  try:
    _add_conformer_to_mol(mol, conformer, force_parse)
    mol = sanitize_mol(mol, sort_alphabetically, remove_hydrogens)
  except (
      ValueError,
      rd_chem.KekulizeException,
      rd_chem.AtomValenceException,
  ) as e:
    raise MolFromMmcifError from e

  return mol


def mol_to_ccd_cif(
    mol: rd_chem.Mol,
    component_id: str,
    pdbx_smiles: str | None = None,
    include_hydrogens: bool = True,
) -> cif_dict.CifDict:
  """Creates a CCD-like mmcif data block from an rdkit Mol object.

  Only a subset of associated mmcif fields is populated, but that is
  sufficient for further usage, e.g. in featurization code.

  Atom names can be specified via `atom_name` property. For atoms with
  unspecified value of that property, the name is assigned based on element type
  and the order in the Mol object.

  If the Mol object has associated conformers, atom positions from the first of
  them will be populated in the resulting mmcif file.

  Args:
     mol: An rdkit molecule.
     component_id: Name of the molecule to use in the resulting mmcif. That is
       equivalent to CCD code.
     pdbx_smiles: If specified, the value will be used to populate
       `_chem_comp.pdbx_smiles`.
     include_hydrogens: Whether to include atom and bond data involving
       hydrogens.

  Returns:
     An mmcif data block corresponding for the given rdkit molecule.

  Raises:
    UnsupportedMolBond: When a molecule contains a bond that can't be
      represented with mmcif.
  """
  mol = rd_chem.Mol(mol)
  if include_hydrogens:
    mol = rd_chem.AddHs(mol)
  rd_chem.Kekulize(mol)

  if mol.GetNumConformers() > 0:
    ideal_conformer = mol.GetConformer(0).GetPositions()
    ideal_conformer = np.vectorize(lambda x: f'{x:.3f}')(ideal_conformer)
  else:
    # No data will be populated in the resulting mmcif if the molecule doesn't
    # have any conformers attached to it.
    ideal_conformer = None

  mol_cif = collections.defaultdict(list)
  mol_cif['data_'] = [component_id]
  mol_cif['_chem_comp.id'] = [component_id]
  if pdbx_smiles:
    mol_cif['_chem_comp.pdbx_smiles'] = [pdbx_smiles]

  mol = assign_atom_names_from_graph(mol, keep_existing_names=True)

  for atom_idx, atom in enumerate(mol.GetAtoms()):
    element = atom.GetSymbol()
    if not include_hydrogens and element in ('H', 'D'):
      continue

    mol_cif['_chem_comp_atom.comp_id'].append(component_id)
    mol_cif['_chem_comp_atom.atom_id'].append(atom.GetProp('atom_name'))
    mol_cif['_chem_comp_atom.type_symbol'].append(atom.GetSymbol().upper())
    mol_cif['_chem_comp_atom.charge'].append(str(atom.GetFormalCharge()))
    if ideal_conformer is not None:
      coords = ideal_conformer[atom_idx]
      mol_cif['_chem_comp_atom.pdbx_model_Cartn_x_ideal'].append(coords[0])
      mol_cif['_chem_comp_atom.pdbx_model_Cartn_y_ideal'].append(coords[1])
      mol_cif['_chem_comp_atom.pdbx_model_Cartn_z_ideal'].append(coords[2])

  for bond in mol.GetBonds():
    atom1 = bond.GetBeginAtom()
    atom2 = bond.GetEndAtom()
    if not include_hydrogens and (
        atom1.GetSymbol() in ('H', 'D') or atom2.GetSymbol() in ('H', 'D')
    ):
      continue
    mol_cif['_chem_comp_bond.comp_id'].append(component_id)
    mol_cif['_chem_comp_bond.atom_id_1'].append(
        bond.GetBeginAtom().GetProp('atom_name')
    )
    mol_cif['_chem_comp_bond.atom_id_2'].append(
        bond.GetEndAtom().GetProp('atom_name')
    )
    try:
      bond_type = bond.GetBondType()
      # Older versions of RDKit did not have a DATIVE bond type. Convert it to
      # SINGLE to match the AF3 training setup.
      if bond_type == rd_chem.BondType.DATIVE:
        bond_type = rd_chem.BondType.SINGLE
      mol_cif['_chem_comp_bond.value_order'].append(
          _RDKIT_BOND_TYPE_TO_MMCIF[bond_type]
      )
      mol_cif['_chem_comp_bond.pdbx_stereo_config'].append(
          _RDKIT_BOND_STEREO_TO_MMCIF[bond.GetStereo()]
      )
    except KeyError as e:
      raise UnsupportedMolBondError from e
    mol_cif['_chem_comp_bond.pdbx_aromatic_flag'].append(
        'Y' if bond.GetIsAromatic() else 'N'
    )

  return cif_dict.CifDict(mol_cif)


def _format_atom_name(atom_name: str, atom_type: str) -> str:
  """Formats an atom name to fit in the four characters specified in PDB.

  See for example the following note on atom name formatting in PDB files:
  https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html#note1

  Args:
    atom_name: The unformatted atom name.
    atom_type: The atom element symbol.

  Returns:
    formatted_atom_name: The formatted 4-character atom name.
  """
  atom_name = atom_name.strip()
  atom_type = atom_type.strip().upper()
  if len(atom_name) == 1:
    return atom_name.rjust(2).ljust(4)
  elif len(atom_name) == 2:
    if atom_name == atom_type:
      return atom_name.ljust(4)
    return atom_name.center(4)
  elif len(atom_name) == 3:
    if atom_name[:2] == atom_type:
      return atom_name.ljust(4)
    return atom_name.rjust(4)
  elif len(atom_name) == 4:
    return atom_name
  else:
    raise ValueError(
        f'Atom name `{atom_name}` has more than four characters '
        'or is an empty string.'
    )


def parse_atom_data(
    mol_cif: cif_dict.CifDict | Mapping[str, Sequence[str]], force_parse: bool
) -> tuple[Sequence[str], Sequence[str], Sequence[int], Sequence[str]]:
  """Parses atoms. If force_parse is True, fix deuterium and missing charge."""
  atom_types = [t.capitalize() for t in mol_cif['_chem_comp_atom.type_symbol']]
  atom_names = mol_cif['_chem_comp_atom.atom_id']
  atom_charges = mol_cif['_chem_comp_atom.charge']
  atom_leaving_flags = ['?'] * len(atom_names)
  if '_chem_comp_atom.pdbx_leaving_atom_flag' in mol_cif:
    atom_leaving_flags = mol_cif['_chem_comp_atom.pdbx_leaving_atom_flag']

  if force_parse:
    # Replace missing charges with 0.
    atom_charges = [charge if charge != '?' else '0' for charge in atom_charges]
    # Deuterium for hydrogen.
    atom_types = [type_ if type_ != 'D' else 'H' for type_ in atom_types]

  atom_charges = [int(atom_charge) for atom_charge in atom_charges]
  return atom_names, atom_types, atom_charges, atom_leaving_flags


def parse_bond_data(
    mol_cif: cif_dict.CifDict | Mapping[str, Sequence[str]], force_parse: bool
) -> tuple[
    Sequence[str], Sequence[str], Sequence[rd_chem.BondType], Sequence[bool]
]:
  """Parses bond data. If force_parse is True, ignore missing aromatic flags."""
  # The bond table isn't present if there are no bonds. Use [] in that case.
  begin_atoms = mol_cif.get('_chem_comp_bond.atom_id_1', [])
  end_atoms = mol_cif.get('_chem_comp_bond.atom_id_2', [])
  orders = mol_cif.get('_chem_comp_bond.value_order', [])
  bond_types = [_RDKIT_MMCIF_TO_BOND_TYPE[order] for order in orders]

  try:
    aromatic_flags = mol_cif.get('_chem_comp_bond.pdbx_aromatic_flag', [])
    is_aromatic = [{'Y': True, 'N': False}[flag] for flag in aromatic_flags]
  except KeyError:
    if force_parse:
      # Set them all to not aromatic.
      is_aromatic = [False for _ in begin_atoms]
    else:
      raise

  return begin_atoms, end_atoms, bond_types, is_aromatic


def _parse_ideal_conformer(mol_cif: cif_dict.CifDict) -> rd_chem.Conformer:
  """Builds a conformer containing the ideal coordinates from the CCD.

  Args:
     mol_cif: An mmcif object representing a molecule.

  Returns:
     An rdkit conformer filled with the ideal positions from the mmcif.

  Raises:
     ValueError: if the positions can't be interpreted.
  """
  atom_x = [
      float(x) for x in mol_cif['_chem_comp_atom.pdbx_model_Cartn_x_ideal']
  ]
  atom_y = [
      float(y) for y in mol_cif['_chem_comp_atom.pdbx_model_Cartn_y_ideal']
  ]
  atom_z = [
      float(z) for z in mol_cif['_chem_comp_atom.pdbx_model_Cartn_z_ideal']
  ]
  atom_positions = zip(atom_x, atom_y, atom_z, strict=True)

  conformer = rd_chem.Conformer(len(atom_x))
  for atom_index, atom_position in enumerate(atom_positions):
    conformer.SetAtomPosition(atom_index, atom_position)

  return conformer


def sort_atoms_by_name(mol: rd_chem.Mol) -> rd_chem.Mol:
  """Sorts the atoms in the molecule by their names."""
  atom_names = {
      atom.GetProp('atom_name'): atom.GetIdx() for atom in mol.GetAtoms()
  }

  # Sort the name, int tuples by the names.
  sorted_atom_names = sorted(atom_names.items())

  # Zip these tuples back together to the sorted indices.
  _, new_order = zip(*sorted_atom_names, strict=True)

  # Reorder the molecule.
  # new_order is effectively an argsort of the names.
  return rd_chem.RenumberAtoms(mol, new_order)


def assign_atom_names_from_graph(
    mol: rd_chem.Mol,
    keep_existing_names: bool = False,
) -> rd_chem.Mol:
  """Assigns atom names from the molecular graph.

  The atom name is stored as an atom property 'atom_name', accessible
  with atom.GetProp('atom_name'). If the property is already specified, and
  keep_existing_names is True we keep the original name.

  We traverse the graph in the order of the rdkit atom index and give each atom
  a name equal to '{ELEMENT_TYPE}{INDEX}'. E.g. C5 is the name for the fifth
  unnamed carbon encountered.

  NOTE: A new mol is returned, the original is not changed in place.

  Args:
    mol: Mol object.
    keep_existing_names: If True, atoms that already have the atom_name property
      will keep their assigned names.

  Returns:
    A new mol, with potentially new 'atom_name' properties.
  """
  mol = rd_chem.Mol(mol)

  specified_atom_names = {
      atom.GetProp('atom_name')
      for atom in mol.GetAtoms()
      if atom.HasProp('atom_name') and keep_existing_names
  }

  element_counts = collections.Counter()
  for atom in mol.GetAtoms():
    if not atom.HasProp('atom_name') or not keep_existing_names:
      element = atom.GetSymbol()
      while True:
        element_counts[element] += 1
        # Standardize names by using uppercase element type, as in CCD. Only
        # effects elements with more than one letter, e.g. 'Cl' becomes 'CL'.
        new_name = f'{element.upper()}{element_counts[element]}'
        if new_name not in specified_atom_names:
          break
      atom.SetProp('atom_name', new_name)

  return mol


def get_random_conformer(
    mol: rd_chem.Mol,
    random_seed: int,
    max_iterations: int | None,
    logging_name: str,
) -> rd_chem.Conformer | None:
  """Stochastic conformer search method using V3 ETK."""
  params = rd_all_chem.ETKDGv3()
  params.randomSeed = random_seed
  if max_iterations is not None:  # Override default value.
    params.maxIterations = max_iterations
  mol_copy = rd_chem.Mol(mol)
  try:
    conformer_id = rd_all_chem.EmbedMolecule(mol_copy, params)
    conformer = mol_copy.GetConformer(conformer_id)
  except ValueError:
    logging.warning('Failed to generate conformer for: %s', logging_name)
    conformer = None
  return conformer
