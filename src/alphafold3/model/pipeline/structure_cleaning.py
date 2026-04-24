# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Prepare PDB structure for training or inference."""

from absl import logging
from alphafold3 import structure
from alphafold3.constants import chemical_component_sets
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.model.atom_layout import atom_layout
from alphafold3.model.pipeline import inter_chain_bonds
import numpy as np


def _get_leaving_atom_mask(
    struc: structure.Structure,
    polymer_ligand_bonds: atom_layout.AtomLayout | None,
    ligand_ligand_bonds: atom_layout.AtomLayout | None,
    chain_id: str,
    chain_type: str,
    res_id: int,
    res_name: str,
) -> np.ndarray:
  """Updates a drop_leaving_atoms mask with new leaving atom locations."""
  bonded_atoms = atom_layout.get_bonded_atoms(
      polymer_ligand_bonds,
      ligand_ligand_bonds,
      res_id,
      chain_id,
  )
  # Connect the amino-acids, i.e. remove OXT, HXT and H2.
  drop_atoms = atom_layout.get_link_drop_atoms(
      res_name=res_name,
      chain_type=chain_type,
      is_start_terminus=False,
      is_end_terminus=False,
      bonded_atoms=bonded_atoms,
      drop_ligand_leaving_atoms=True,
  )
  # Default mask where everything is false, which equates to being kept.
  drop_atom_filter_atoms = struc.chain_id != struc.chain_id
  for drop_atom in drop_atoms:
    drop_atom_filter_atom = np.logical_and(
        np.logical_and(
            struc.atom_name == drop_atom,
            struc.chain_id == chain_id,
        ),
        struc.res_id == res_id,
    )
    drop_atom_filter_atoms = np.logical_or(
        drop_atom_filter_atoms, drop_atom_filter_atom
    )
  return drop_atom_filter_atoms


def clean_structure(
    struc: structure.Structure,
    ccd: chemical_components.Ccd,
    *,
    drop_missing_sequence: bool,
    drop_non_standard_atoms: bool,
    filter_waters: bool,
    filter_hydrogens: bool,
    filter_leaving_atoms: bool,
    only_glycan_ligands_for_leaving_atoms: bool,
    covalent_bonds_only: bool,
    remove_polymer_polymer_bonds: bool,
    remove_bad_bonds: bool,
) -> structure.Structure:
  """Returns a cleaned version of the input structure.

  Args:
    struc: Structure to clean.
    ccd: The chemical components dictionary.
    drop_missing_sequence: Whether to drop chains without specified sequences.
    drop_non_standard_atoms: Whether to drop non CCD standard atoms.
    filter_waters: Whether to drop water chains.
    filter_hydrogens: Whether to drop hyrdogen atoms.
    filter_leaving_atoms: Whether to drop leaving atoms based on heuristics.
    only_glycan_ligands_for_leaving_atoms: Whether to only include glycan
      ligands when filtering leaving atoms.
    covalent_bonds_only: Only include covalent bonds.
    remove_polymer_polymer_bonds: Remove polymer-polymer bonds.
    remove_bad_bonds: Whether to remove badly bonded ligands.
  """

  # Drop chains without specified sequences.
  if drop_missing_sequence:
    chains_with_unk_sequence = struc.find_chains_with_unknown_sequence()
    if chains_with_unk_sequence:
      struc = struc.filter_out(chain_id=chains_with_unk_sequence)

  # Drop non-standard atoms
  if drop_non_standard_atoms:
    struc = struc.drop_non_standard_atoms(
        ccd=ccd, drop_unk=False, drop_non_ccd=False
    )

  # Sort chains in "reverse-spreadsheet" order.
  struc = struc.with_sorted_chains

  if filter_hydrogens:
    struc = struc.without_hydrogen()

  if filter_waters:
    struc = struc.filter_out(chain_type=mmcif_names.WATER)

  if filter_leaving_atoms:
    drop_leaving_atoms_all = struc.chain_id != struc.chain_id
    polymer_ligand_bonds = inter_chain_bonds.get_polymer_ligand_bonds(
        struc,
        only_glycan_ligands=only_glycan_ligands_for_leaving_atoms,
    )
    ligand_ligand_bonds = inter_chain_bonds.get_ligand_ligand_bonds(
        struc,
        only_glycan_ligands=only_glycan_ligands_for_leaving_atoms,
    )
    all_glycans = {
        *chemical_component_sets.GLYCAN_OTHER_LIGANDS,
        *chemical_component_sets.GLYCAN_LINKING_LIGANDS,
    }
    # If only glycan ligands and no O1 atoms, we can do parallel drop. If there are no glycan bonds at all, keep O1.
    has_any_glycan_bond = bool(
        ligand_ligand_bonds.atom_name.size or polymer_ligand_bonds.atom_name.size
    )
    if (
        only_glycan_ligands_for_leaving_atoms
        and has_any_glycan_bond
        and (not (ligand_ligand_bonds.atom_name == 'O1').any())
        and (not (polymer_ligand_bonds.atom_name == 'O1').any())
    ):
      drop_leaving_atoms_all = np.logical_and(
          np.isin(struc.atom_name, 'O1'),
          np.isin(struc.res_name, list(all_glycans)),
      )
    else:
      substruct = struc.group_by_residue
      glycan_mask = np.isin(substruct.res_name, list(all_glycans))
      substruct = substruct.filter(glycan_mask)
      # We need to iterate over all glycan residues for this.
      for res in substruct.iter_residues():
        # Only need to do drop leaving atoms for glycans depending on bonds.
        if (res_name := res['res_name']) in all_glycans:
          drop_atom_filter = _get_leaving_atom_mask(
              struc=struc,
              polymer_ligand_bonds=polymer_ligand_bonds,
              ligand_ligand_bonds=ligand_ligand_bonds,
              chain_id=res['chain_id'],
              chain_type=res['chain_type'],
              res_id=res['res_id'],
              res_name=res_name,
          )
          drop_leaving_atoms_all = np.logical_or(
              drop_leaving_atoms_all, drop_atom_filter
          )

    num_atoms_before = struc.num_atoms
    struc = struc.filter_out(drop_leaving_atoms_all)
    num_atoms_after = struc.num_atoms

    if num_atoms_before > num_atoms_after:
      logging.error(
          'Dropped %s atoms from GT struc: chain_id %s res_id %s res_name %s',
          num_atoms_before - num_atoms_after,
          struc.chain_id,
          struc.res_id,
          struc.res_name,
      )

  # Can filter by bond type without having to iterate over bonds.
  if struc.bonds and covalent_bonds_only:
    is_covalent = np.isin(struc.bonds.type, ['covale'])
    if sum(is_covalent) > 0:
      new_bonds = struc.bonds[is_covalent]
    else:
      new_bonds = structure.Bonds.make_empty()
    struc = struc.copy_and_update(bonds=new_bonds)

  # Other bond filters require iterating over individual bonds.
  if struc.bonds and (remove_bad_bonds or remove_polymer_polymer_bonds):
    include_bond = []
    num_pp_bonds = 0
    num_bad_bonds = 0
    for bond in struc.iter_bonds():
      dest_atom = bond.dest_atom
      from_atom = bond.from_atom
      if remove_polymer_polymer_bonds:
        if (
            from_atom['chain_type'] in mmcif_names.POLYMER_CHAIN_TYPES
            and dest_atom['chain_type'] in mmcif_names.POLYMER_CHAIN_TYPES
        ):
          num_pp_bonds += 1
          include_bond.append(False)
          continue
      if remove_bad_bonds:
        dest_coords = np.array(
            [dest_atom['atom_x'], dest_atom['atom_y'], dest_atom['atom_z']]
        )
        from_coords = np.array(
            [from_atom['atom_x'], from_atom['atom_y'], from_atom['atom_z']]
        )
        squared_dist = np.sum(np.square(dest_coords - from_coords))
        squared_threshold = 2.4 * 2.4
        if squared_dist > squared_threshold:
          num_bad_bonds += 1
          include_bond.append(False)
          continue
      include_bond.append(True)
    if sum(include_bond) < len(struc.bonds):
      logging.info(
          'Reducing number of bonds for %s from %s to %s, of which %s are'
          ' polymer-polymer bonds and %s are bad bonds.',
          struc.name,
          len(struc.bonds),
          sum(include_bond),
          num_pp_bonds,
          num_bad_bonds,
      )
      if sum(include_bond) > 0:
        # Need to index bonds with bond keys or arrays of bools with same length
        # as num bonds. In this case, we use array of bools (as elsewhere in the
        # cleaning code).
        new_bonds = struc.bonds[np.array(include_bond, dtype=bool)]
      else:
        new_bonds = structure.Bonds.make_empty()
      struc = struc.copy_and_update(bonds=new_bonds)

  return struc


def create_empty_output_struc_and_layout(
    struc: structure.Structure,
    ccd: chemical_components.Ccd,
    *,
    with_hydrogens: bool = False,
    skip_unk: bool = False,
    polymer_ligand_bonds: atom_layout.AtomLayout | None = None,
    ligand_ligand_bonds: atom_layout.AtomLayout | None = None,
    drop_ligand_leaving_atoms: bool = False,
) -> tuple[structure.Structure, atom_layout.AtomLayout]:
  """Make zero-coordinate structure from all physical residues.

  Args:
    struc: Structure object.
    ccd: The chemical components dictionary.
    with_hydrogens: Whether to keep hydrogen atoms in structure.
    skip_unk: Whether to remove unknown residues from structure.
    polymer_ligand_bonds: Bond information for polymer-ligand pairs.
    ligand_ligand_bonds: Bond information for ligand-ligand pairs.
    drop_ligand_leaving_atoms: Flag for handling leaving atoms for ligands.

  Returns:
    Tuple of structure with all bonds, physical residues and coordinates set to
    0 and a flat atom layout of empty structure.
  """
  bonded_atom_pairs = []
  if polymer_ligand_bonds:
    for chain_ids, res_ids, atom_names in zip(
        polymer_ligand_bonds.chain_id,
        polymer_ligand_bonds.res_id,
        polymer_ligand_bonds.atom_name,
        strict=True,
    ):
      bonded_atom_pairs.append((
          (chain_ids[0], res_ids[0], atom_names[0]),
          (chain_ids[1], res_ids[1], atom_names[1]),
      ))
  if ligand_ligand_bonds:
    for chain_ids, res_ids, atom_names in zip(
        ligand_ligand_bonds.chain_id,
        ligand_ligand_bonds.res_id,
        ligand_ligand_bonds.atom_name,
        strict=True,
    ):
      bonded_atom_pairs.append((
          (chain_ids[0], res_ids[0], atom_names[0]),
          (chain_ids[1], res_ids[1], atom_names[1]),
      ))
  residues = atom_layout.residues_from_structure(
      struc, include_missing_residues=True
  )

  flat_output_layout = atom_layout.make_flat_atom_layout(
      residues,
      ccd=ccd,
      with_hydrogens=with_hydrogens,
      skip_unk_residues=skip_unk,
      polymer_ligand_bonds=polymer_ligand_bonds,
      ligand_ligand_bonds=ligand_ligand_bonds,
      drop_ligand_leaving_atoms=drop_ligand_leaving_atoms,
  )

  empty_output_struc = atom_layout.make_structure(
      flat_layout=flat_output_layout,
      atom_coords=np.zeros((flat_output_layout.shape[0], 3)),
      name=struc.name,
      atom_b_factors=None,
      all_physical_residues=residues,
  )
  if bonded_atom_pairs:
    empty_output_struc = empty_output_struc.add_bonds(
        bonded_atom_pairs, bond_type=mmcif_names.COVALENT_BOND
    )

  return empty_output_struc, flat_output_layout
