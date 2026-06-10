# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Functions for handling inter-chain bonds."""

from collections.abc import Collection
import functools
from typing import Final, NamedTuple
from alphafold3 import structure
from alphafold3.constants import chemical_component_sets
from alphafold3.constants import mmcif_names
from alphafold3.model.atom_layout import atom_layout
import numpy as np


BOND_THRESHOLD_GLYCANS_ANGSTROM: Final[float] = 1.7
# See https://pubs.acs.org/doi/10.1021/ja010331r for P-P atom bond distances.
BOND_THRESHOLD_ALL_ANGSTROM: Final[float] = 2.4


class BondAtomArrays(NamedTuple):
  chain_id: np.ndarray
  chain_type: np.ndarray
  res_id: np.ndarray
  res_name: np.ndarray
  atom_name: np.ndarray
  coords: np.ndarray


def _get_bond_atom_arrays(
    struc: structure.Structure, bond_atom_indices: np.ndarray
) -> BondAtomArrays:
  return BondAtomArrays(
      chain_id=struc.chain_id[bond_atom_indices],
      chain_type=struc.chain_type[bond_atom_indices],
      res_id=struc.res_id[bond_atom_indices],
      res_name=struc.res_name[bond_atom_indices],
      atom_name=struc.atom_name[bond_atom_indices],
      coords=struc.coords[..., bond_atom_indices, :],
  )


@functools.lru_cache(maxsize=1)
def get_polymer_ligand_and_ligand_ligand_bonds(
    struct: structure.Structure,
    only_glycan_ligands: bool,
    allow_multiple_bonds_per_atom: bool,
) -> tuple[atom_layout.AtomLayout, atom_layout.AtomLayout]:
  """Return polymer-ligand & ligand-ligand inter-residue bonds.

  Args:
    struct: Structure object to extract bonds from.
    only_glycan_ligands: Whether to only include glycans in ligand category.
    allow_multiple_bonds_per_atom: If not allowed, we greedily choose the first
      bond seen per atom and discard the remaining on each atom..

  Returns:
    polymer_ligand, ligand_ligand_bonds: Each object is an AtomLayout object
    [num_bonds, 2] for the bond-defining atoms.
  """
  if only_glycan_ligands:
    allowed_res_names = list({
        *chemical_component_sets.GLYCAN_OTHER_LIGANDS,
        *chemical_component_sets.GLYCAN_LINKING_LIGANDS,
    })
  else:
    allowed_res_names = None
  all_bonds = get_bond_layout(
      bond_threshold=BOND_THRESHOLD_GLYCANS_ANGSTROM
      if only_glycan_ligands
      else BOND_THRESHOLD_ALL_ANGSTROM,
      struct=struct,
      allowed_chain_types1=list({
          *mmcif_names.LIGAND_CHAIN_TYPES,
          *mmcif_names.POLYMER_CHAIN_TYPES,
      }),
      allowed_chain_types2=list(mmcif_names.LIGAND_CHAIN_TYPES),
      allowed_res_names=allowed_res_names,
      allow_multiple_bonds_per_atom=allow_multiple_bonds_per_atom,
  )
  ligand_ligand_bonds_mask = np.isin(
      all_bonds.chain_type, list(mmcif_names.LIGAND_CHAIN_TYPES)
  )
  polymer_ligand_bonds_mask = np.isin(
      all_bonds.chain_type, list(mmcif_names.POLYMER_CHAIN_TYPES)
  )
  polymer_ligand_bonds_mask = np.logical_and(
      ligand_ligand_bonds_mask.any(axis=1),
      polymer_ligand_bonds_mask.any(axis=1),
  )
  ligand_ligand_bonds = all_bonds[ligand_ligand_bonds_mask.all(axis=1)]
  polymer_ligand_bonds = all_bonds[polymer_ligand_bonds_mask]
  return polymer_ligand_bonds, ligand_ligand_bonds


def _remove_multi_bonds(
    bond_layout: atom_layout.AtomLayout,
) -> atom_layout.AtomLayout:
  """Remove instances greedily."""
  uids = {}
  keep_indx = []
  for chain_id, res_id, atom_name in zip(
      bond_layout.chain_id,
      bond_layout.res_id,
      bond_layout.atom_name,
      strict=True,
  ):
    key1 = (chain_id[0], res_id[0], atom_name[0])
    key2 = (chain_id[1], res_id[1], atom_name[1])
    keep_indx.append(bool(key1 not in uids) and bool(key2 not in uids))
    if key1 not in uids:
      uids[key1] = None
    if key2 not in uids:
      uids[key2] = None
  return bond_layout[np.array(keep_indx, dtype=bool)]


@functools.lru_cache(maxsize=1)
def get_ligand_ligand_bonds(
    struct: structure.Structure,
    only_glycan_ligands: bool,
    allow_multiple_bonds_per_atom: bool = False,
) -> atom_layout.AtomLayout:
  """Return ligand-ligand inter-residue bonds.

  Args:
    struct: Structure object to extract bonds from.
    only_glycan_ligands: Whether to only include glycans in ligand category.
    allow_multiple_bonds_per_atom: If not allowed, we greedily choose the first
      bond seen per atom and discard the remaining on each atom.

  Returns:
    bond_layout: AtomLayout object [num_bonds, 2] for the bond-defining atoms.
  """
  if only_glycan_ligands:
    allowed_res_names = list({
        *chemical_component_sets.GLYCAN_OTHER_LIGANDS,
        *chemical_component_sets.GLYCAN_LINKING_LIGANDS,
    })
  else:
    allowed_res_names = None
  return get_bond_layout(
      bond_threshold=BOND_THRESHOLD_GLYCANS_ANGSTROM
      if only_glycan_ligands
      else BOND_THRESHOLD_ALL_ANGSTROM,
      struct=struct,
      allowed_chain_types1=list(mmcif_names.LIGAND_CHAIN_TYPES),
      allowed_chain_types2=list(mmcif_names.LIGAND_CHAIN_TYPES),
      allowed_res_names=allowed_res_names,
      allow_multiple_bonds_per_atom=allow_multiple_bonds_per_atom,
  )


@functools.lru_cache(maxsize=1)
def get_polymer_ligand_bonds(
    struct: structure.Structure,
    only_glycan_ligands: bool,
    allow_multiple_bonds_per_atom: bool = False,
    bond_threshold: float | None = None,
) -> atom_layout.AtomLayout:
  """Return polymer-ligand interchain bonds.

  Args:
    struct: Structure object to extract bonds from.
    only_glycan_ligands: Whether to only include glycans in ligand category.
    allow_multiple_bonds_per_atom: If not allowed, we greedily choose the first
      bond seen per atom and discard the remaining on each atom.
    bond_threshold: Euclidean distance of max allowed bond.

  Returns:
    bond_layout: AtomLayout object [num_bonds, 2] for the bond-defining atoms.
  """
  if only_glycan_ligands:
    allowed_res_names = list({
        *chemical_component_sets.GLYCAN_OTHER_LIGANDS,
        *chemical_component_sets.GLYCAN_LINKING_LIGANDS,
    })
  else:
    allowed_res_names = None
  if bond_threshold is None:
    if only_glycan_ligands:
      bond_threshold = BOND_THRESHOLD_GLYCANS_ANGSTROM
    else:
      bond_threshold = BOND_THRESHOLD_ALL_ANGSTROM
  return get_bond_layout(
      bond_threshold=bond_threshold,
      struct=struct,
      allowed_chain_types1=list(mmcif_names.POLYMER_CHAIN_TYPES),
      allowed_chain_types2=list(mmcif_names.LIGAND_CHAIN_TYPES),
      allowed_res_names=allowed_res_names,
      allow_multiple_bonds_per_atom=allow_multiple_bonds_per_atom,
  )


def get_bond_layout(
    bond_threshold: float = BOND_THRESHOLD_ALL_ANGSTROM,
    *,
    struct: structure.Structure,
    allowed_chain_types1: Collection[str],
    allowed_chain_types2: Collection[str],
    include_bond_types: Collection[str] = ('covale',),
    allowed_res_names: Collection[str] | None = None,
    allow_multiple_bonds_per_atom: bool,
) -> atom_layout.AtomLayout:
  """Get bond_layout for all bonds between two sets of chain types.

  There is a mask (all_mask) that runs through this script, and each bond pair
  needs to maintain a True across all conditions in order to be preserved at the
  end, otherwise the bond pair has invalidated a condition with a False and is
  removed entirely. Note, we remove oxygen atom bonds as they are an edge case
  that causes issues with scoring, due to multiple waters bonding with single
  residues.

  Args:
    bond_threshold: Maximum bond distance in Angstrom.
    struct: Structure object to extract bonds from.
    allowed_chain_types1: One end of the bonds must be an atom with one of these
      chain types.
    allowed_chain_types2: The other end of the bond must be an atom with one of
      these chain types.
    include_bond_types: Only include bonds with specified type e.g. hydrog,
      metalc, covale, disulf.
    allowed_res_names: Further restricts from chain_types. Either end of the
      bonds must be an atom part of these res_names. If none all will be
      accepted after chain and bond type filtering.
    allow_multiple_bonds_per_atom: If not allowed, we greedily choose the first
      bond seen per atom and discard the remaining on each atom.

  Returns:
    bond_layout: AtomLayout object [num_bonds, 2] for the bond-defining atoms.
  """
  if not struct.bonds:
    return atom_layout.AtomLayout(
        atom_name=np.empty((0, 2), dtype=object),
        res_id=np.empty((0, 2), dtype=int),
        res_name=np.empty((0, 2), dtype=object),
        chain_id=np.empty((0, 2), dtype=object),
        chain_type=np.empty((0, 2), dtype=object),
        atom_element=np.empty((0, 2), dtype=object),
    )
  from_atom_idxs, dest_atom_idxs = struct.bonds.get_atom_indices(
      struct.atom_key
  )
  from_atoms = _get_bond_atom_arrays(struct, from_atom_idxs)
  dest_atoms = _get_bond_atom_arrays(struct, dest_atom_idxs)
  # Chain type
  chain_mask = np.logical_or(
      np.logical_and(
          np.isin(
              from_atoms.chain_type,
              allowed_chain_types1,
          ),
          np.isin(
              dest_atoms.chain_type,
              allowed_chain_types2,
          ),
      ),
      np.logical_and(
          np.isin(
              from_atoms.chain_type,
              allowed_chain_types2,
          ),
          np.isin(
              dest_atoms.chain_type,
              allowed_chain_types1,
          ),
      ),
  )
  if allowed_res_names:
    # Res type
    res_mask = np.logical_or(
        np.isin(from_atoms.res_name, allowed_res_names),
        np.isin(dest_atoms.res_name, allowed_res_names),
    )
    # All mask
    all_mask = np.logical_and(chain_mask, res_mask)
  else:
    all_mask = chain_mask
  # Bond type mask
  type_mask = np.isin(struct.bonds.type, list(include_bond_types))
  np.logical_and(all_mask, type_mask, out=all_mask)
  # Bond length check. Work in square length to avoid taking many square roots.
  bond_length_squared = np.square(from_atoms.coords - dest_atoms.coords).sum(
      axis=1
  )
  bond_threshold_squared = bond_threshold * bond_threshold
  np.logical_and(
      all_mask, bond_length_squared < bond_threshold_squared, out=all_mask
  )
  # Inter-chain and inter-residue bonds for ligands
  ligand_types = list(mmcif_names.LIGAND_CHAIN_TYPES)
  is_ligand = np.logical_or(
      np.isin(
          from_atoms.chain_type,
          ligand_types,
      ),
      np.isin(
          dest_atoms.chain_type,
          ligand_types,
      ),
  )
  res_id_differs = from_atoms.res_id != dest_atoms.res_id
  chain_id_differs = from_atoms.chain_id != dest_atoms.chain_id
  is_inter_res = np.logical_or(res_id_differs, chain_id_differs)
  is_inter_ligand_res = np.logical_and(is_inter_res, is_ligand)
  is_inter_chain_not_ligand = np.logical_and(chain_id_differs, ~is_ligand)
  # If ligand then inter-res & inter-chain bonds, otherwise inter-chain only.
  combined_allowed_bonds = np.logical_or(
      is_inter_chain_not_ligand, is_inter_ligand_res
  )
  np.logical_and(all_mask, combined_allowed_bonds, out=all_mask)
  bond_layout = atom_layout.AtomLayout(
      atom_name=np.stack(
          [
              from_atoms.atom_name[all_mask],
              dest_atoms.atom_name[all_mask],
          ],
          axis=1,
          dtype=object,
      ),
      res_id=np.stack(
          [from_atoms.res_id[all_mask], dest_atoms.res_id[all_mask]],
          axis=1,
          dtype=int,
      ),
      chain_id=np.stack(
          [
              from_atoms.chain_id[all_mask],
              dest_atoms.chain_id[all_mask],
          ],
          axis=1,
          dtype=object,
      ),
  )
  if not allow_multiple_bonds_per_atom:
    bond_layout = _remove_multi_bonds(bond_layout)
  return atom_layout.fill_in_optional_fields(
      bond_layout,
      reference_atoms=atom_layout.atom_layout_from_structure(struct),
  )
