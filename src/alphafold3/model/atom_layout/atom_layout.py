# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Helper functions for different atom layouts and conversion between them."""

import collections
from collections.abc import Mapping, Sequence
import dataclasses
import types
from typing import Any, TypeAlias

from alphafold3 import structure
from alphafold3.constants import atom_types
from alphafold3.constants import chemical_component_sets
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.constants import residue_names
from alphafold3.data.tools import rdkit_utils
from alphafold3.structure import chemical_components as struc_chem_comps
import jax.numpy as jnp
import numpy as np
from rdkit import Chem


xnp_ndarray: TypeAlias = np.ndarray | jnp.ndarray  # pylint: disable=invalid-name
NumpyIndex: TypeAlias = Any


@dataclasses.dataclass(frozen=True)
class AtomLayout:
  """Atom layout in a fixed shape (usually 1-dim or 2-dim).

  Examples for atom layouts are atom37, atom14, and similar.
  All members are np.ndarrays with the same shape, e.g.
  - [num_atoms]
  - [num_residues, max_atoms_per_residue]
  - [num_fragments, max_fragments_per_residue]
  All string arrays should have dtype=object to avoid pitfalls with Numpy's
  fixed-size strings

  Attributes:
    atom_name: np.ndarray of str: atom names (e.g. 'CA', 'NE2'), padding
      elements have an empty string (''), None or any other value, that maps to
      False for .astype(bool). mmCIF field: _atom_site.label_atom_id.
    res_id: np.ndarray of int: residue index (usually starting from 1) padding
      elements can have an arbitrary value. mmCIF field:
      _atom_site.label_seq_id.
    chain_id: np.ndarray of str: chain names (e.g. 'A', 'B') padding elements
      can have an arbitrary value. mmCIF field: _atom_site.label_seq_id.
    atom_element: np.ndarray of str: atom elements (e.g. 'C', 'N', 'O'), padding
      elements have an empty string (''), None or any other value, that maps to
      False for .astype(bool). mmCIF field: _atom_site.type_symbol.
    res_name: np.ndarray of str: residue names (e.g. 'ARG', 'TRP') padding
      elements can have an arbitrary value. mmCIF field:
      _atom_site.label_comp_id.
    chain_type: np.ndarray of str: chain types (e.g. 'polypeptide(L)'). padding
      elements can have an arbitrary value. mmCIF field: _entity_poly.type OR
      _entity.type (for non-polymers).
    shape: shape of the layout (just returns atom_name.shape)
  """

  atom_name: np.ndarray
  res_id: np.ndarray
  chain_id: np.ndarray
  atom_element: np.ndarray | None = None
  res_name: np.ndarray | None = None
  chain_type: np.ndarray | None = None

  def __post_init__(self):
    """Assert all arrays have the same shape."""
    attribute_names = (
        'atom_name',
        'atom_element',
        'res_name',
        'res_id',
        'chain_id',
        'chain_type',
    )
    _assert_all_arrays_have_same_shape(
        obj=self,
        expected_shape=self.atom_name.shape,
        attribute_names=attribute_names,
    )
    # atom_name must have dtype object, such that we can convert it to bool to
    # obtain the mask
    if self.atom_name.dtype != object:
      raise ValueError(
          'atom_name must have dtype object, such that it can '
          'be converted converted to bool to obtain the mask'
      )

  def __getitem__(self, key: NumpyIndex) -> 'AtomLayout':
    return AtomLayout(
        atom_name=self.atom_name[key],
        res_id=self.res_id[key],
        chain_id=self.chain_id[key],
        atom_element=(
            self.atom_element[key] if self.atom_element is not None else None
        ),
        res_name=(self.res_name[key] if self.res_name is not None else None),
        chain_type=(
            self.chain_type[key] if self.chain_type is not None else None
        ),
    )

  def __eq__(self, other: 'AtomLayout') -> bool:
    if not np.array_equal(self.atom_name, other.atom_name):
      return False

    mask = self.atom_name.astype(bool)
    # Check essential fields.
    for field in ('res_id', 'chain_id'):
      my_arr = getattr(self, field)
      other_arr = getattr(other, field)
      if not np.array_equal(my_arr[mask], other_arr[mask]):
        return False

    # Check optional fields.
    for field in ('atom_element', 'res_name', 'chain_type'):
      my_arr = getattr(self, field)
      other_arr = getattr(other, field)
      if (
          my_arr is not None
          and other_arr is not None
          and not np.array_equal(my_arr[mask], other_arr[mask])
      ):
        return False

    return True

  def copy_and_pad_to(self, shape: tuple[int, ...]) -> 'AtomLayout':
    """Copies and pads the layout to the requested shape.

    Args:
      shape: new shape for the atom layout

    Returns:
      a copy of the atom layout padded to the requested shape

    Raises:
      ValueError: incompatible shapes.
    """
    if len(shape) != len(self.atom_name.shape):
      raise ValueError(
          f'Incompatible shape {shape}. Current layout has shape {self.shape}.'
      )
    if any(new < old for old, new in zip(self.atom_name.shape, shape)):
      raise ValueError(
          "Can't pad to a smaller shape. Current layout has shape "
          f'{self.shape} and you requested shape {shape}.'
      )
    pad_width = [
        (0, new - old) for old, new in zip(self.atom_name.shape, shape)
    ]
    pad_val = np.array('', dtype=object)
    return AtomLayout(
        atom_name=np.pad(self.atom_name, pad_width, constant_values=pad_val),
        res_id=np.pad(self.res_id, pad_width, constant_values=0),
        chain_id=np.pad(self.chain_id, pad_width, constant_values=pad_val),
        atom_element=(
            np.pad(self.atom_element, pad_width, constant_values=pad_val)
            if self.atom_element is not None
            else None
        ),
        res_name=(
            np.pad(self.res_name, pad_width, constant_values=pad_val)
            if self.res_name is not None
            else None
        ),
        chain_type=(
            np.pad(self.chain_type, pad_width, constant_values=pad_val)
            if self.chain_type is not None
            else None
        ),
    )

  def to_array(self) -> np.ndarray:
    """Stacks the fields to a numpy array with shape (6, <layout_shape>).

    Creates a pure numpy array of type `object` by stacking the 6 fields of the
    AtomLayout, i.e. (atom_name, atom_element, res_name, res_id, chain_id,
    chain_type). This method together with from_array() provides an easy way to
    apply pure numpy methods like np.concatenate() to `AtomLayout`s.

    Returns:
      np.ndarray of object with shape (6, <layout_shape>), e.g.
      array([['N', 'CA', 'C', ..., 'CB', 'CG', 'CD'],
       ['N', 'C', 'C', ..., 'C', 'C', 'C'],
       ['LEU', 'LEU', 'LEU', ..., 'PRO', 'PRO', 'PRO'],
       [1, 1, 1, ..., 403, 403, 403],
       ['A', 'A', 'A', ..., 'D', 'D', 'D'],
       ['polypeptide(L)', 'polypeptide(L)', ..., 'polypeptide(L)']],
      dtype=object)
    """
    if (
        self.atom_element is None
        or self.res_name is None
        or self.chain_type is None
    ):
      raise ValueError('All optional fields need to be present.')

    return np.stack(dataclasses.astuple(self), axis=0)

  @classmethod
  def from_array(cls, arr: np.ndarray) -> 'AtomLayout':
    """Creates an AtomLayout object from a numpy array with shape (6, ...).

    see also to_array()
    Args:
      arr: np.ndarray of object with shape (6, <layout_shape>)

    Returns:
      AtomLayout object with shape (<layout_shape>)
    """
    if arr.shape[0] != 6:
      raise ValueError(
          'Given array must have shape (6, ...) to match the 6 fields of '
          'AtomLayout (atom_name, atom_element, res_name, res_id, chain_id, '
          f'chain_type). Your array has {arr.shape=}'
      )
    return cls(*arr)

  @property
  def shape(self) -> tuple[int, ...]:
    return self.atom_name.shape


@dataclasses.dataclass(frozen=True)
class Residues:
  """List of residues with meta data.

  Attributes:
    res_name: np.ndarray of str [num_res], e.g. 'ARG', 'TRP'
    res_id: np.ndarray of int [num_res]
    chain_id: np.ndarray of str [num_res], e.g. 'A', 'B'
    chain_type: np.ndarray of str [num_res], e.g. 'polypeptide(L)'
    is_start_terminus: np.ndarray of bool [num_res]
    is_end_terminus: np.ndarray of bool [num_res]
    deprotonation: (optional) np.ndarray of set() [num_res], e.g. {'HD1', 'HE2'}
    smiles_string: (optional) np.ndarray of str [num_res], e.g. 'Cc1ccccc1'
    shape: shape of the layout (just returns res_name.shape)
  """

  res_name: np.ndarray
  res_id: np.ndarray
  chain_id: np.ndarray
  chain_type: np.ndarray
  is_start_terminus: np.ndarray
  is_end_terminus: np.ndarray
  deprotonation: np.ndarray | None = None
  smiles_string: np.ndarray | None = None

  def __post_init__(self):
    """Assert all arrays are 1D have the same shape."""
    attribute_names = (
        'res_name',
        'res_id',
        'chain_id',
        'chain_type',
        'is_start_terminus',
        'is_end_terminus',
        'deprotonation',
        'smiles_string',
    )
    _assert_all_arrays_have_same_shape(
        obj=self,
        expected_shape=(self.res_name.shape[0],),
        attribute_names=attribute_names,
    )

  def __getitem__(self, key: NumpyIndex) -> 'Residues':
    return Residues(
        res_name=self.res_name[key],
        res_id=self.res_id[key],
        chain_id=self.chain_id[key],
        chain_type=self.chain_type[key],
        is_start_terminus=self.is_start_terminus[key],
        is_end_terminus=self.is_end_terminus[key],
        deprotonation=(
            self.deprotonation[key] if self.deprotonation is not None else None
        ),
        smiles_string=(
            self.smiles_string[key] if self.smiles_string is not None else None
        ),
    )

  def __eq__(self, other: 'Residues') -> bool:
    return all(
        np.array_equal(getattr(self, field.name), getattr(other, field.name))
        for field in dataclasses.fields(self)
    )

  @property
  def shape(self) -> tuple[int, ...]:
    return self.res_name.shape


@dataclasses.dataclass(frozen=True)
class GatherInfo:
  """Gather indices to translate from one atom layout to another.

  All members are np or jnp ndarray (usually 1-dim or 2-dim) with the same
  shape, e.g.
  - [num_atoms]
  - [num_residues, max_atoms_per_residue]
  - [num_fragments, max_fragments_per_residue]

  Attributes:
    gather_idxs: np or jnp ndarray of int: gather indices into a flattened array
    gather_mask: np or jnp ndarray of bool: mask for resulting array
    input_shape: np or jnp ndarray of int: the shape of the unflattened input
      array
    shape: output shape. Just returns gather_idxs.shape
  """

  gather_idxs: xnp_ndarray
  gather_mask: xnp_ndarray
  input_shape: xnp_ndarray

  def __post_init__(self):
    if self.gather_mask.shape != self.gather_idxs.shape:
      raise ValueError(
          'All arrays must have the same shape. Got\n'
          f'gather_idxs.shape = {self.gather_idxs.shape}\n'
          f'gather_mask.shape = {self.gather_mask.shape}\n'
      )

  def __getitem__(self, key: NumpyIndex) -> 'GatherInfo':
    return GatherInfo(
        gather_idxs=self.gather_idxs[key],
        gather_mask=self.gather_mask[key],
        input_shape=self.input_shape,
    )

  @property
  def shape(self) -> tuple[int, ...]:
    return self.gather_idxs.shape

  def as_np_or_jnp(self, xnp: types.ModuleType) -> 'GatherInfo':
    return GatherInfo(
        gather_idxs=xnp.array(self.gather_idxs),
        gather_mask=xnp.array(self.gather_mask),
        input_shape=xnp.array(self.input_shape),
    )

  def as_dict(
      self,
      key_prefix: str | None = None,
  ) -> dict[str, xnp_ndarray]:
    prefix = f'{key_prefix}:' if key_prefix else ''
    return {
        prefix + 'gather_idxs': self.gather_idxs,
        prefix + 'gather_mask': self.gather_mask,
        prefix + 'input_shape': self.input_shape,
    }

  @classmethod
  def from_dict(
      cls,
      d: Mapping[str, xnp_ndarray],
      key_prefix: str | None = None,
  ) -> 'GatherInfo':
    """Creates GatherInfo from a given dictionary."""
    prefix = f'{key_prefix}:' if key_prefix else ''
    return cls(
        gather_idxs=d[prefix + 'gather_idxs'],
        gather_mask=d[prefix + 'gather_mask'],
        input_shape=d[prefix + 'input_shape'],
    )


def fill_in_optional_fields(
    minimal_atom_layout: AtomLayout,
    reference_atoms: AtomLayout,
) -> AtomLayout:
  """Fill in the optional fields (atom_element, res_name, chain_type).

  Extracts the optional fields (atom_element, res_name, chain_type) from a
  flat reference layout and fills them into the fields from this layout.

  Args:
    minimal_atom_layout: An AtomLayout that only contains the essential fields
      (atom_name, res_id, chain_id).
    reference_atoms: A flat layout that contains all fields for all atoms.

  Returns:
    An AtomLayout that contains all fields.

  Raises:
    ValueError: Reference atoms layout is not flat.
    ValueError: Missing atoms in reference.
  """
  if len(reference_atoms.shape) > 1:
    raise ValueError('Only flat layouts are supported as reference.')
  ref_to_self = compute_gather_idxs(
      source_layout=reference_atoms, target_layout=minimal_atom_layout
  )
  atom_mask = minimal_atom_layout.atom_name.astype(bool)
  missing_atoms_mask = atom_mask & ~ref_to_self.gather_mask
  if np.any(missing_atoms_mask):
    raise ValueError(
        f'{np.sum(missing_atoms_mask)} missing atoms in reference: '
        f'{minimal_atom_layout[missing_atoms_mask]}'
    )

  def _convert_str_array(gather: GatherInfo, arr: np.ndarray):
    output = arr[gather.gather_idxs]
    output[~gather.gather_mask] = ''
    return output

  return dataclasses.replace(
      minimal_atom_layout,
      atom_element=_convert_str_array(
          ref_to_self, reference_atoms.atom_element
      ),
      res_name=_convert_str_array(ref_to_self, reference_atoms.res_name),
      chain_type=_convert_str_array(ref_to_self, reference_atoms.chain_type),
  )


def guess_deprotonation(residues: Residues) -> Residues:
  """Convenience function to create a plausible deprotonation field.

  Assumes a pH of 7 and always prefers HE2 over HD1 for HIS.
  Args:
    residues: a Residues object without a depronotation field

  Returns:
    a Residues object with a depronotation field
  """
  num_residues = residues.res_name.shape[0]
  deprotonation = np.empty(num_residues, dtype=object)
  deprotonation_at_ph7 = {
      'ASP': 'HD2',
      'GLU': 'HE2',
      'HIS': 'HD1',
  }
  for idx, res_name in enumerate(residues.res_name):
    deprotonation[idx] = set()
    if res_name in deprotonation_at_ph7:
      deprotonation[idx].add(deprotonation_at_ph7[res_name])
    if residues.is_end_terminus[idx]:
      deprotonation[idx].add('HXT')

  return dataclasses.replace(residues, deprotonation=deprotonation)


def atom_layout_from_structure(
    struct: structure.Structure,
    *,
    fix_non_standard_polymer_res: bool = False,
) -> AtomLayout:
  """Extract AtomLayout from a Structure."""

  if not fix_non_standard_polymer_res:
    return AtomLayout(
        atom_name=np.array(struct.atom_name, dtype=object),
        atom_element=np.array(struct.atom_element, dtype=object),
        res_name=np.array(struct.res_name, dtype=object),
        res_id=np.array(struct.res_id, dtype=int),
        chain_id=np.array(struct.chain_id, dtype=object),
        chain_type=np.array(struct.chain_type, dtype=object),
    )

  # Target lists.
  target_atom_names = []
  target_atom_elements = []
  target_res_ids = []
  target_res_names = []
  target_chain_ids = []
  target_chain_types = []

  for atom in struct.iter_atoms():
    target_atom_names.append(atom['atom_name'])
    target_atom_elements.append(atom['atom_element'])
    target_res_ids.append(atom['res_id'])
    target_chain_ids.append(atom['chain_id'])
    target_chain_types.append(atom['chain_type'])
    if mmcif_names.is_standard_polymer_type(atom['chain_type']):
      fixed_res_name = mmcif_names.fix_non_standard_polymer_res(
          res_name=atom['res_name'], chain_type=atom['chain_type']
      )
      target_res_names.append(fixed_res_name)
    else:
      target_res_names.append(atom['res_name'])

  return AtomLayout(
      atom_name=np.array(target_atom_names, dtype=object),
      atom_element=np.array(target_atom_elements, dtype=object),
      res_name=np.array(target_res_names, dtype=object),
      res_id=np.array(target_res_ids, dtype=int),
      chain_id=np.array(target_chain_ids, dtype=object),
      chain_type=np.array(target_chain_types, dtype=object),
  )


def residues_from_structure(
    struct: structure.Structure,
    *,
    include_missing_residues: bool = True,
    fix_non_standard_polymer_res: bool = False,
) -> Residues:
  """Create a Residues object from a Structure object."""

  def _get_smiles(res_name):
    """Get SMILES string from chemical components."""
    smiles = None
    if (
        struct.chemical_components_data is not None
        and struct.chemical_components_data.chem_comp is not None
        and struct.chemical_components_data.chem_comp.get(res_name)
    ):
      smiles = struct.chemical_components_data.chem_comp[res_name].pdbx_smiles
    return smiles

  res_names_per_chain = struct.chain_res_name_sequence(
      include_missing_residues=include_missing_residues,
      fix_non_standard_polymer_res=fix_non_standard_polymer_res,
  )
  res_name = []
  res_id = []
  chain_id = []
  chain_type = []
  smiles = []
  is_start_terminus = []
  for c in struct.iter_chains():
    if include_missing_residues:
      this_res_ids = [id for (_, id) in struct.all_residues[c['chain_id']]]
    else:
      this_res_ids = [
          r['res_id']
          for r in struct.iter_residues()
          if r['chain_id'] == c['chain_id']
      ]
    fixed_res_names = res_names_per_chain[c['chain_id']]
    assert len(this_res_ids) == len(
        fixed_res_names
    ), f'{len(this_res_ids)} != {len(fixed_res_names)}'
    this_start_res_id = min(min(this_res_ids), 1)
    this_is_start_terminus = [r == this_start_res_id for r in this_res_ids]
    smiles.extend([_get_smiles(res_name) for res_name in fixed_res_names])
    num_res = len(fixed_res_names)
    res_name.extend(fixed_res_names)
    res_id.extend(this_res_ids)
    chain_id.extend([c['chain_id']] * num_res)
    chain_type.extend([c['chain_type']] * num_res)
    is_start_terminus.extend(this_is_start_terminus)
  res_name = np.array(res_name, dtype=object)
  res_id = np.array(res_id, dtype=int)
  chain_id = np.array(chain_id, dtype=object)
  chain_type = np.array(chain_type, dtype=object)
  smiles = np.array(smiles, dtype=object)
  is_start_terminus = np.array(is_start_terminus, dtype=bool)

  res_uid_to_idx = {
      uid: idx for idx, uid in enumerate(zip(chain_id, res_id, strict=True))
  }

  # Start terminus indicates whether residue index is 1 and chain is polymer.
  is_polymer = np.isin(chain_type, tuple(mmcif_names.POLYMER_CHAIN_TYPES))
  is_start_terminus = is_start_terminus & is_polymer

  # Start also indicates whether amino acid is attached to H2 or proline to H.
  start_terminus_atom_index = np.nonzero(
      (struct.chain_type == mmcif_names.PROTEIN_CHAIN)
      & (
          (struct.atom_name == 'H2')
          | ((struct.atom_name == 'H') & (struct.res_name == 'PRO'))
      )
  )[0]

  # Translate atom idx to residue idx to assign start terminus.
  for atom_idx in start_terminus_atom_index:
    res_uid = (struct.chain_id[atom_idx], struct.res_id[atom_idx])
    res_idx = res_uid_to_idx[res_uid]
    is_start_terminus[res_idx] = True

  # Infer end terminus: Check for OXT, or in case of
  # include_missing_residues==True for the last residue of the chain.
  num_all_residues = res_name.shape[0]
  is_end_terminus = np.zeros(num_all_residues, dtype=bool)
  end_term_atom_idxs = np.nonzero(struct.atom_name == 'OXT')[0]
  for atom_idx in end_term_atom_idxs:
    res_uid = (struct.chain_id[atom_idx], struct.res_id[atom_idx])
    res_idx = res_uid_to_idx[res_uid]
    is_end_terminus[res_idx] = True

  if include_missing_residues:
    for idx in range(num_all_residues - 1):
      if is_polymer[idx] and chain_id[idx] != chain_id[idx + 1]:
        is_end_terminus[idx] = True
    if (num_all_residues > 0) and is_polymer[-1]:
      is_end_terminus[-1] = True

  # Infer (de-)protonation: Only if hydrogens are given.
  num_hydrogens = np.sum(
      (struct.atom_element == 'H') & (struct.chain_type == 'polypeptide(L)')
  )
  if num_hydrogens > 0:
    deprotonation = np.empty(num_all_residues, dtype=object)
    all_atom_uids = set(
        zip(struct.chain_id, struct.res_id, struct.atom_name, strict=True)
    )
    for idx in range(num_all_residues):
      deprotonation[idx] = set()
      check_hydrogens = set()
      if is_end_terminus[idx]:
        check_hydrogens.add('HXT')
      if res_name[idx] in atom_types.PROTONATION_HYDROGENS:
        check_hydrogens.update(atom_types.PROTONATION_HYDROGENS[res_name[idx]])
      for hydrogen in check_hydrogens:
        if (chain_id[idx], res_id[idx], hydrogen) not in all_atom_uids:
          deprotonation[idx].add(hydrogen)
  else:
    deprotonation = None

  return Residues(
      res_name=res_name,
      res_id=res_id,
      chain_id=chain_id,
      chain_type=chain_type,
      is_start_terminus=is_start_terminus.astype(bool),
      is_end_terminus=is_end_terminus,
      deprotonation=deprotonation,
      smiles_string=smiles,
  )


def get_link_drop_atoms(
    res_name: str,
    chain_type: str,
    *,
    is_start_terminus: bool,
    is_end_terminus: bool,
    bonded_atoms: set[str],
    drop_ligand_leaving_atoms: bool = False,
) -> set[str]:
  """Returns set of atoms that are dropped when this res_name gets linked.

  Args:
    res_name: residue name, e.g. 'ARG'
    chain_type: chain_type, e.g. 'polypeptide(L)'
    is_start_terminus: whether the residue is the n-terminus
    is_end_terminus: whether the residue is the c-terminus
    bonded_atoms: Names of atoms coming off this residue.
    drop_ligand_leaving_atoms: Flag to switch on/off leaving atoms for ligands.

  Returns:
    Set of atoms that are dropped when this amino acid gets linked.
  """
  drop_atoms = set()
  if chain_type == mmcif_names.PROTEIN_CHAIN:
    if res_name == 'PRO':
      if not is_start_terminus:
        drop_atoms.update({'H', 'H2', 'H3'})
      if not is_end_terminus:
        drop_atoms.update({'OXT', 'HXT'})
    else:
      if not is_start_terminus:
        drop_atoms.update({'H2', 'H3'})
      if not is_end_terminus:
        drop_atoms.update({'OXT', 'HXT'})
  elif chain_type in mmcif_names.NUCLEIC_ACID_CHAIN_TYPES:
    if not is_start_terminus:
      drop_atoms.update({'OP3'})
  elif (
      drop_ligand_leaving_atoms and chain_type in mmcif_names.LIGAND_CHAIN_TYPES
  ):
    if res_name in {
        *chemical_component_sets.GLYCAN_OTHER_LIGANDS,
        *chemical_component_sets.GLYCAN_LINKING_LIGANDS,
    }:
      if bonded_atoms and 'O1' not in bonded_atoms:
        drop_atoms.update({'O1'})
  return drop_atoms


def get_bonded_atoms(
    polymer_ligand_bonds: AtomLayout,
    ligand_ligand_bonds: AtomLayout,
    res_id: int,
    chain_id: str,
) -> set[str]:
  """Finds the res_name on the opposite end of the bond, if a bond exists.

  Args:
    polymer_ligand_bonds: Bond information for polymer-ligand pairs.
    ligand_ligand_bonds: Bond information for ligand-ligand pairs.
    res_id: residue id in question.
    chain_id: chain id of residue in question.

  Returns:
    res_name of bonded atom.
  """
  bonded_atoms = set()
  if polymer_ligand_bonds:
    # Filter before searching to speed this up.
    bond_idx = np.logical_and(
        polymer_ligand_bonds.res_id == res_id,
        polymer_ligand_bonds.chain_id == chain_id,
    ).any(axis=1)
    relevant_polymer_bonds = polymer_ligand_bonds[bond_idx]
    for atom_names, res_ids, chain_ids in zip(
        relevant_polymer_bonds.atom_name,
        relevant_polymer_bonds.res_id,
        relevant_polymer_bonds.chain_id,
    ):
      if (res_ids[0], chain_ids[0]) == (res_id, chain_id):
        bonded_atoms.add(atom_names[0])
      elif (res_ids[1], chain_ids[1]) == (res_id, chain_id):
        bonded_atoms.add(atom_names[1])
  if ligand_ligand_bonds:
    bond_idx = np.logical_and(
        ligand_ligand_bonds.res_id == res_id,
        ligand_ligand_bonds.chain_id == chain_id,
    ).any(axis=1)
    relevant_ligand_bonds = ligand_ligand_bonds[bond_idx]
    for atom_names, res_ids, chain_ids in zip(
        relevant_ligand_bonds.atom_name,
        relevant_ligand_bonds.res_id,
        relevant_ligand_bonds.chain_id,
    ):
      if (res_ids[0], chain_ids[0]) == (res_id, chain_id):
        bonded_atoms.add(atom_names[0])
      elif (res_ids[1], chain_ids[1]) == (res_id, chain_id):
        bonded_atoms.add(atom_names[1])
  return bonded_atoms


def make_flat_atom_layout(
    residues: Residues,
    ccd: chemical_components.Ccd,
    polymer_ligand_bonds: AtomLayout | None = None,
    ligand_ligand_bonds: AtomLayout | None = None,
    *,
    with_hydrogens: bool = False,
    skip_unk_residues: bool = True,
    drop_ligand_leaving_atoms: bool = False,
) -> AtomLayout:
  """Make a flat atom layout for given residues.

  Create a flat layout from a `Residues` object. The required atoms for each
  amino acid type are taken from the CCD, hydrogens and oxygens are dropped to
  make the linked residues. Terminal OXT's and protonation state for the
  hydrogens come from the `Residues` object.

  Args:
    residues: a `Residues` object.
    ccd: The chemical components dictionary.
    polymer_ligand_bonds: Bond information for polymer-ligand pairs.
    ligand_ligand_bonds: Bond information for ligand-ligand pairs.
    with_hydrogens: whether to create hydrogens
    skip_unk_residues: whether to skip 'UNK' resides -- default is True to be
      compatible with the rest of AlphaFold that does not predict atoms for
      unknown residues
    drop_ligand_leaving_atoms: Flag to switch on/ off leaving atoms for ligands.

  Returns:
    an `AtomLayout` object
  """
  num_res = residues.res_name.shape[0]

  # Target lists.
  target_atom_names = []
  target_atom_elements = []
  target_res_ids = []
  target_res_names = []
  target_chain_ids = []
  target_chain_types = []

  for idx in range(num_res):
    # skip 'UNK' residues if requested
    if (
        skip_unk_residues
        and residues.res_name[idx] in residue_names.UNKNOWN_TYPES
    ):
      continue

    # Get the atoms for this residue type from CCD.
    if ccd.get(residues.res_name[idx]):
      res_atoms = struc_chem_comps.get_all_atoms_in_entry(
          ccd=ccd, res_name=residues.res_name[idx]
      )
      atom_names_elements = list(
          zip(
              res_atoms['_chem_comp_atom.atom_id'],
              res_atoms['_chem_comp_atom.type_symbol'],
              strict=True,
          )
      )
    elif residues.smiles_string[idx]:
      # Get atoms from RDKit via SMILES.
      mol = Chem.MolFromSmiles(residues.smiles_string[idx])
      if mol is None:
        raise ValueError(
            f'Failed to construct RDKit Mol for {residues.res_name[idx]} from'
            f' SMILES string: {residues.smiles_string[idx]} . This is likely'
            ' due to an issue with the SMILES string. Note that the userCCD'
            ' input format provides an alternative way to define custom'
            ' molecules directly without RDKit or SMILES.'
        )
      mol = rdkit_utils.assign_atom_names_from_graph(mol)
      atom_names_elements = [
          (a.GetProp('atom_name'), a.GetSymbol()) for a in mol.GetAtoms()
      ]
    else:
      raise ValueError(
          f'{residues.res_name[idx]} not found in CCD and no SMILES string'
      )

    # Remove hydrogens if requested.
    if not with_hydrogens:
      atom_names_elements = [
          (n, e) for n, e in atom_names_elements if (e != 'H' and e != 'D')
      ]
    bonded_atoms = get_bonded_atoms(
        polymer_ligand_bonds,
        ligand_ligand_bonds,
        residues.res_id[idx],
        residues.chain_id[idx],
    )
    # Connect the amino-acids, i.e. remove OXT, HXT and H2.
    drop_atoms = get_link_drop_atoms(
        res_name=residues.res_name[idx],
        chain_type=residues.chain_type[idx],
        is_start_terminus=residues.is_start_terminus[idx],
        is_end_terminus=residues.is_end_terminus[idx],
        bonded_atoms=bonded_atoms,
        drop_ligand_leaving_atoms=drop_ligand_leaving_atoms,
    )

    # If deprotonation info is available, remove the specific atoms.
    if residues.deprotonation is not None:
      drop_atoms.update(residues.deprotonation[idx])

    atom_names_elements = [
        (n, e) for n, e in atom_names_elements if n not in drop_atoms
    ]

    # Append the found atoms to the target lists.
    target_atom_names.extend([n for n, _ in atom_names_elements])
    target_atom_elements.extend([e for _, e in atom_names_elements])
    num_atoms = len(atom_names_elements)
    target_res_names.extend([residues.res_name[idx]] * num_atoms)
    target_res_ids.extend([residues.res_id[idx]] * num_atoms)
    target_chain_ids.extend([residues.chain_id[idx]] * num_atoms)
    target_chain_types.extend([residues.chain_type[idx]] * num_atoms)

  return AtomLayout(
      atom_name=np.array(target_atom_names, dtype=object),
      atom_element=np.array(target_atom_elements, dtype=object),
      res_name=np.array(target_res_names, dtype=object),
      res_id=np.array(target_res_ids, dtype=int),
      chain_id=np.array(target_chain_ids, dtype=object),
      chain_type=np.array(target_chain_types, dtype=object),
  )


def compute_gather_idxs(
    *,
    source_layout: AtomLayout,
    target_layout: AtomLayout,
    fill_value: int = 0,
) -> GatherInfo:
  """Produce gather indices and mask to convert from source layout to target."""
  source_uid_to_idx = {
      uid: idx
      for idx, uid in enumerate(
          zip(
              source_layout.chain_id.ravel(),
              source_layout.res_id.ravel(),
              source_layout.atom_name.ravel(),
              strict=True,
          )
      )
  }
  gather_idxs = []
  gather_mask = []
  for uid in zip(
      target_layout.chain_id.ravel(),
      target_layout.res_id.ravel(),
      target_layout.atom_name.ravel(),
      strict=True,
  ):
    if uid in source_uid_to_idx:
      gather_idxs.append(source_uid_to_idx[uid])
      gather_mask.append(True)
    else:
      gather_idxs.append(fill_value)
      gather_mask.append(False)
  target_shape = target_layout.atom_name.shape
  return GatherInfo(
      gather_idxs=np.array(gather_idxs, dtype=int).reshape(target_shape),
      gather_mask=np.array(gather_mask, dtype=bool).reshape(target_shape),
      input_shape=np.array(source_layout.atom_name.shape),
  )


def convert(
    gather_info: GatherInfo,
    arr: xnp_ndarray,
    *,
    layout_axes: tuple[int, ...] = (0,),
) -> xnp_ndarray:
  """Convert an array from one atom layout to another."""
  # Translate negative indices to the corresponding positives.
  layout_axes = tuple(i if i >= 0 else i + arr.ndim for i in layout_axes)

  # Ensure that layout_axes are continuous.
  layout_axes_begin = layout_axes[0]
  layout_axes_end = layout_axes[-1] + 1

  if layout_axes != tuple(range(layout_axes_begin, layout_axes_end)):
    raise ValueError(f'layout_axes must be continuous. Got {layout_axes}.')
  layout_shape = arr.shape[layout_axes_begin:layout_axes_end]

  # Ensure that the layout shape is compatible
  # with the gather_info. I.e. the first axis size must be equal or greater
  # than the gather_info.input_shape, and all subsequent axes sizes must match.
  if (len(layout_shape) != gather_info.input_shape.size) or (
      isinstance(gather_info.input_shape, np.ndarray)
      and (
          (layout_shape[0] < gather_info.input_shape[0])
          or (np.any(layout_shape[1:] != gather_info.input_shape[1:]))
      )
  ):
    raise ValueError(
        'Input array layout axes are incompatible. You specified layout '
        f'axes {layout_axes} with an input array of shape {arr.shape}, but '
        f'the gather info expects shape {gather_info.input_shape}. '
        'Your first axis size must be equal or greater than the '
        'gather_info.input_shape, and all subsequent axes sizes must '
        'match.'
    )

  # Compute the shape of the input array with flattened layout.
  batch_shape = arr.shape[:layout_axes_begin]
  features_shape = arr.shape[layout_axes_end:]
  arr_flattened_shape = batch_shape + (np.prod(layout_shape),) + features_shape

  # Flatten input array and perform the gather.
  arr_flattened = arr.reshape(arr_flattened_shape)
  if layout_axes_begin == 0:
    out_arr = arr_flattened[gather_info.gather_idxs, ...]
  elif layout_axes_begin == 1:
    out_arr = arr_flattened[:, gather_info.gather_idxs, ...]
  elif layout_axes_begin == 2:
    out_arr = arr_flattened[:, :, gather_info.gather_idxs, ...]
  elif layout_axes_begin == 3:
    out_arr = arr_flattened[:, :, :, gather_info.gather_idxs, ...]
  elif layout_axes_begin == 4:
    out_arr = arr_flattened[:, :, :, :, gather_info.gather_idxs, ...]
  else:
    raise ValueError(
        'Only 4 batch axes supported. If you need more, the code '
        'is easy to extend.'
    )

  # Broadcast the mask and apply it.
  broadcasted_mask_shape = (
      (1,) * len(batch_shape)
      + gather_info.gather_mask.shape
      + (1,) * len(features_shape)
  )
  out_arr *= gather_info.gather_mask.reshape(broadcasted_mask_shape)
  return out_arr


def make_structure(
    flat_layout: AtomLayout,
    atom_coords: np.ndarray,
    name: str,
    *,
    atom_b_factors: np.ndarray | None = None,
    all_physical_residues: Residues | None = None,
) -> structure.Structure:
  """Returns a Structure from a flat layout and atom coordinates.

  The provided flat_layout must be 1-dim and must not contain any padding
  elements. The flat_layout.atom_name must conform to the OpenMM/CCD standard
  and must not contain deuterium.

  Args:
    flat_layout: flat 1-dim AtomLayout without pading elements
    atom_coords: np.ndarray of float, shape (num_atoms, 3)
    name: str: the name (usually PDB id), e.g. '1uao'
    atom_b_factors: np.ndarray of float, shape (num_atoms,) or None. If None,
      they will be set to all zeros.
    all_physical_residues: a Residues object that contains all physically
      existing residues, i.e. also those residues that have no resolved atoms.
      This is common in experimental structures, but also appears in predicted
      structures for 'UNK' or other non-standard residue types, where the model
      does not predict coordinates. This will be used to create the
      `all_residues` field of the structure object.
  """

  if flat_layout.atom_name.ndim != 1 or not np.all(
      flat_layout.atom_name.astype(bool)
  ):
    raise ValueError(
        'flat_layout must be 1-dim and must not contain anypadding element'
    )
  if (
      flat_layout.atom_element is None
      or flat_layout.res_name is None
      or flat_layout.chain_type is None
  ):
    raise ValueError('All optional fields must be present.')

  if atom_b_factors is None:
    atom_b_factors = np.zeros(atom_coords.shape[:-1])

  if all_physical_residues is not None:
    # Create the all_residues field from a Residues object
    # (unfortunately there is no central place to keep the chain_types in
    # the structure class, so we drop it here)
    all_residues = collections.defaultdict(list)
    for chain_id, res_id, res_name in zip(
        all_physical_residues.chain_id,
        all_physical_residues.res_id,
        all_physical_residues.res_name,
        strict=True,
    ):
      all_residues[chain_id].append((res_name, res_id))
  else:
    # Create the all_residues field from the flat_layout
    all_residues = collections.defaultdict(list)
    if flat_layout.chain_id.shape[0] > 0:
      all_residues[flat_layout.chain_id[0]].append(
          (flat_layout.res_name[0], flat_layout.res_id[0])
      )
      for i in range(1, flat_layout.shape[0]):
        if (
            flat_layout.chain_id[i] != flat_layout.chain_id[i - 1]
            or flat_layout.res_name[i] != flat_layout.res_name[i - 1]
            or flat_layout.res_id[i] != flat_layout.res_id[i - 1]
        ):
          all_residues[flat_layout.chain_id[i]].append(
              (flat_layout.res_name[i], flat_layout.res_id[i])
          )

  return structure.from_atom_arrays(
      name=name,
      all_residues=dict(all_residues),
      chain_id=flat_layout.chain_id,
      chain_type=flat_layout.chain_type,
      res_id=flat_layout.res_id.astype(np.int32),
      res_name=flat_layout.res_name,
      atom_name=flat_layout.atom_name,
      atom_element=flat_layout.atom_element,
      atom_x=atom_coords[..., 0],
      atom_y=atom_coords[..., 1],
      atom_z=atom_coords[..., 2],
      atom_b_factor=atom_b_factors,
  )


def _assert_all_arrays_have_same_shape(
    *,
    obj: AtomLayout | Residues | GatherInfo,
    expected_shape: tuple[int, ...],
    attribute_names: Sequence[str],
) -> None:
  """Checks that given attributes of the object have the expected shape."""
  attribute_shapes_description = []
  all_shapes_are_valid = True

  for attribute_name in attribute_names:
    attribute = getattr(obj, attribute_name)

    if attribute is None:
      attribute_shape = None
    else:
      attribute_shape = attribute.shape

    if attribute_shape is not None and expected_shape != attribute_shape:
      all_shapes_are_valid = False

    attribute_shape_name = attribute_name + '.shape'
    attribute_shapes_description.append(
        f'{attribute_shape_name:25} = {attribute_shape}'
    )

  if not all_shapes_are_valid:
    raise ValueError(
        f'All arrays must have the same shape ({expected_shape=}). Got\n'
        + '\n'.join(attribute_shapes_description)
    )
