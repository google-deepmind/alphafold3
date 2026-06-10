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

"""Module for parsing various data sources and producing Structures."""

from collections.abc import Collection, Mapping, MutableMapping, Sequence
import dataclasses
import datetime
import enum
import functools
import itertools
from typing import TypeAlias

from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.constants import residue_names
from alphafold3.cpp import mmcif_utils
from alphafold3.cpp import string_array
from alphafold3.structure import bioassemblies
from alphafold3.structure import bonds
from alphafold3.structure import chemical_components as struc_chem_comps
from alphafold3.structure import mmcif
from alphafold3.structure import structure
from alphafold3.structure import structure_tables
import numpy as np


ChainIndex: TypeAlias = int
ResIndex: TypeAlias = int
AtomName: TypeAlias = str
BondAtomId: TypeAlias = tuple[ChainIndex, ResIndex, AtomName]

_INSERTION_CODE_REMAP: Mapping[str, str] = {'.': '?'}


class NoAtomsError(Exception):
  """Raise when the mmCIF does not have any atoms."""


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class BondIndices:
  from_indices: list[int]
  dest_indices: list[int]


@enum.unique
class ModelID(enum.Enum):
  """Values for specifying model IDs when parsing."""

  FIRST = 1  # The first model in the file.
  ALL = 2  # All models in the file.


@enum.unique
class SequenceFormat(enum.Enum):
  """The possible formats for an input sequence."""

  FASTA = 'fasta'  # One-letter code used in FASTA.
  CCD_CODES = 'ccd_codes'  # Multiple-letter chemical components dictionary ids.
  LIGAND_SMILES = 'ligand_smiles'  # SMILES string defining a molecule.


def _create_bond_lookup(
    bonded_atom_pairs: Sequence[tuple[BondAtomId, BondAtomId]],
) -> Mapping[tuple[ChainIndex, ResIndex], Mapping[AtomName, BondIndices]]:
  """Creates maps to help find bonds during a loop over residues."""
  bond_lookup = {}
  for bond_i, (from_atom_id, dest_atom_id) in enumerate(bonded_atom_pairs):
    from_chain_i, from_res_i, from_atom_name = from_atom_id
    dest_chain_i, dest_res_i, dest_atom_name = dest_atom_id
    bonds_by_from_atom_name = bond_lookup.setdefault(
        (from_chain_i, from_res_i), {}
    )
    bonds_by_dest_atom_name = bond_lookup.setdefault(
        (dest_chain_i, dest_res_i), {}
    )
    bonds_by_from_atom_name.setdefault(
        from_atom_name, BondIndices(from_indices=[], dest_indices=[])
    ).from_indices.append(bond_i)
    bonds_by_dest_atom_name.setdefault(
        dest_atom_name, BondIndices(from_indices=[], dest_indices=[])
    ).dest_indices.append(bond_i)
  return bond_lookup


def _get_atom_element(
    ccd: chemical_components.Ccd, res_name: str, atom_name: str
) -> str:
  type_symbol = chemical_components.type_symbol(
      ccd, res_name=res_name, atom_name=atom_name
  )
  return type_symbol or '?'


def _get_representative_atom(
    ccd: chemical_components.Ccd,
    res_name: str,
    chain_type: str,
    sequence_format: SequenceFormat,
) -> tuple[str, str]:
  match sequence_format:
    case SequenceFormat.CCD_CODES:
      atom_name = _get_first_non_leaving_atom(ccd=ccd, res_name=res_name)
      atom_element = _get_atom_element(
          ccd=ccd, res_name=res_name, atom_name=atom_name
      )
      return atom_name, atom_element
    case SequenceFormat.LIGAND_SMILES:
      return '', '?'
    case SequenceFormat.FASTA:
      if chain_type in mmcif_names.PEPTIDE_CHAIN_TYPES:
        return 'CA', 'C'
      if chain_type in mmcif_names.NUCLEIC_ACID_CHAIN_TYPES:
        return "C1'", 'C'
      else:
        raise ValueError(chain_type)
    case _:
      raise ValueError(sequence_format)


@functools.lru_cache(maxsize=128)
def _get_first_non_leaving_atom(
    ccd: chemical_components.Ccd, res_name: str
) -> str:
  """Returns first definitely non-leaving atom if exists, as a stand-in."""
  all_atoms = struc_chem_comps.get_all_atoms_in_entry(ccd, res_name=res_name)[
      '_chem_comp_atom.atom_id'
  ]
  representative_atom = all_atoms[0]
  if representative_atom == 'O1' and len(all_atoms) > 1:
    representative_atom = all_atoms[1]
  return representative_atom


def _add_ligand_to_chem_comp(
    chem_comp: MutableMapping[str, struc_chem_comps.ChemCompEntry],
    ligand_id: str,
    ligand_smiles: str,
):
  """Adds a ligand to chemical components. Raises ValueError on mismatch."""
  new_entry = struc_chem_comps.ChemCompEntry(
      type='non-polymer', pdbx_smiles=ligand_smiles
  )

  existing_entry = chem_comp.get(ligand_id)
  if existing_entry is None:
    chem_comp[ligand_id] = new_entry
  elif existing_entry != new_entry:
    raise ValueError(
        f'Mismatching data for ligand {ligand_id}: '
        f'{new_entry} != {existing_entry}'
    )


def _get_first_model_id(cif: mmcif.Mmcif) -> str:
  """Returns cheaply the first model ID from the mmCIF."""
  return cif.get_array(
      '_atom_site.pdbx_PDB_model_num', dtype=object, gather=slice(1)
  )[0]


def _get_str_model_id(
    cif: mmcif.Mmcif,
    model_id: ModelID | int,
) -> str:
  """Converts a user-specified model_id argument into a string."""
  match model_id:
    case int():
      str_model_id = str(model_id)
    case enum.Enum():
      # We compare the enum's value attribute since regular enum comparison
      # breaks when adhoc importing.
      match model_id.value:
        case ModelID.FIRST.value:
          try:
            str_model_id = _get_first_model_id(cif)
          except IndexError as e:
            raise NoAtomsError(
                'The mmCIF does not have any atoms or'
                ' _atom_site.pdbx_PDB_model_num is missing.'
            ) from e
        case ModelID.ALL.value:
          str_model_id = ''
        case _:
          raise ValueError(
              f'Model ID {model_id} with value {model_id.value} not recognized.'
          )
    case _:
      raise ValueError(
          f'Model ID {model_id} with type {type(model_id)} not recognized.'
      )
  return str_model_id


def _parse_bonds(
    cif: mmcif.Mmcif,
    atom_key: np.ndarray,
    model_id: str,
) -> bonds.Bonds:
  """Returns the bonds table extracted from the mmCIF.

  Args:
    cif: The raw mmCIF to extract the bond information from.
    atom_key: A numpy array defining atom key for each atom in _atom_site. Note
      that the atom key must be computed before resolving alt-locs since this
      function operates on the raw mmCIF!
    model_id: The ID of the model to get bonds for.
  """
  if '_struct_conn.id' not in cif:
    # This is the category key item for the _struct_conn table, therefore
    # we use it to determine whether to parse bond info.
    return bonds.Bonds.make_empty()
  from_atom, dest_atom = mmcif.get_bond_atom_indices(cif, model_id)
  from_atom = np.array(from_atom, dtype=np.int64)
  dest_atom = np.array(dest_atom, dtype=np.int64)
  num_bonds = from_atom.shape[0]
  bond_key = np.arange(num_bonds, dtype=np.int64)
  bond_type = cif.get_array('_struct_conn.conn_type_id', dtype=object)
  if '_struct_conn.pdbx_role' in cif:  # This column isn't always present.
    bond_role = cif.get_array('_struct_conn.pdbx_role', dtype=object)
  else:
    bond_role = np.full((num_bonds,), '?', dtype=object)

  bonds_mask = np.ones((num_bonds,), dtype=bool)
  # Symmetries other than 1_555 imply the atom is not part of the asymmetric
  # unit, and therefore this is a bond that only exists in the expanded
  # bioassembly.
  # We do not currently support parsing these types of bonds.
  if '_struct_conn.ptnr1_symmetry' in cif:
    ptnr1_symmetry = cif.get_array('_struct_conn.ptnr1_symmetry', dtype=object)
    np.logical_and(bonds_mask, ptnr1_symmetry == '1_555', out=bonds_mask)
  if '_struct_conn.ptnr2_symmetry' in cif:
    ptnr2_symmetry = cif.get_array('_struct_conn.ptnr2_symmetry', dtype=object)
    np.logical_and(bonds_mask, ptnr2_symmetry == '1_555', out=bonds_mask)
  # Remove bonds that involve atoms that are not part of the structure,
  # e.g. waters if include_water=False. In a rare case this also removes invalid
  # bonds that are indicated by a key that is set to _atom_site size.
  np.logical_and(bonds_mask, np.isin(from_atom, atom_key), out=bonds_mask)
  np.logical_and(bonds_mask, np.isin(dest_atom, atom_key), out=bonds_mask)
  return bonds.Bonds(
      key=bond_key[bonds_mask],
      type=bond_type[bonds_mask],
      role=bond_role[bonds_mask],
      from_atom_key=from_atom[bonds_mask],
      dest_atom_key=dest_atom[bonds_mask],
  )


@dataclasses.dataclass(frozen=True, slots=True)
class _MmcifHeader:
  name: str
  resolution: float | None
  release_date: datetime.date | None
  structure_method: str | None
  bioassembly_data: bioassemblies.BioassemblyData | None
  chemical_components_data: struc_chem_comps.ChemicalComponentsData | None


def _get_mmcif_header(
    cif: mmcif.Mmcif,
    fix_mse: bool,
    fix_unknown_dna: bool,
) -> _MmcifHeader:
  """Extract header fields from an mmCIF object."""
  entry_id = cif.get('_entry.id')
  name = entry_id[0] if entry_id else cif.get_data_name()
  resolution = mmcif.get_resolution(cif)

  release_date = mmcif.get_release_date(cif)
  if release_date is not None:
    release_date = datetime.date.fromisoformat(release_date)

  experiments = cif.get('_exptl.method')
  structure_method = ','.join(experiments) if experiments else None

  try:
    bioassembly_data = bioassemblies.BioassemblyData.from_mmcif(cif)
  except bioassemblies.MissingBioassemblyDataError:
    bioassembly_data = None

  try:
    chemical_components_data = (
        struc_chem_comps.ChemicalComponentsData.from_mmcif(
            cif, fix_mse=fix_mse, fix_unknown_dna=fix_unknown_dna
        )
    )
  except struc_chem_comps.MissingChemicalComponentsDataError:
    chemical_components_data = None

  return _MmcifHeader(
      name=name,
      resolution=resolution,
      release_date=release_date,
      structure_method=structure_method,
      bioassembly_data=bioassembly_data,
      chemical_components_data=chemical_components_data,
  )


def from_parsed_mmcif(
    mmcif_object: mmcif.Mmcif,
    *,
    name: str | None = None,
    fix_mse_residues: bool = False,
    fix_arginines: bool = False,
    fix_unknown_dna: bool = False,
    include_water: bool = False,
    include_other: bool = False,
    include_bonds: bool = False,
    model_id: int | ModelID = ModelID.FIRST,
) -> structure.Structure:
  """Construct a Structure from a parsed mmCIF object.

  This function is called by `from_mmcif` but can be useful when an mmCIF has
  already been parsed e.g. to extract extra information from the header before
  then converting to Structure for further manipulation.

  Args:
    mmcif_object: A parsed mmcif.Mmcif object.
    name: Optional name for the structure. If not provided, the name will be
      taken from the mmCIF data_ field.
    fix_mse_residues: If True, selenium atom sites (SE) in selenomethionine
      (MSE) residues will be changed to sulphur atom sites (SD). This is because
      methionine (MET) residues are often replaced with MSE to aid X-Ray
      crystallography. If False, the SE MSE atom sites won't be modified.
    fix_arginines: If True, NH1 and NH2 in arginine will be swapped if needed so
      that NH1 is always closer to CD than NH2. If False, no atom sites in
      arginine will be touched. Note that HH11, HH12, HH21, HH22 are fixed too.
    fix_unknown_dna: If True, residues with name N in DNA chains will have their
      res_name replaced with DN. Atoms are not changed.
    include_water: If True, water (HOH) molecules will be parsed. Water
      molecules may be grouped into chains, where number of residues > 1. Water
      molecules are usually grouped into chains but do not necessarily all share
      the same chain ID.
    include_other: If True, all other atoms that are not included by any of the
      above parameters will be included. This covers e.g. "polypeptide(D)" and
      "macrolide" entities, as well as all other non-standard types.
    include_bonds: If True, bond information will be parsed from the mmCIF and
      stored in the Structure.
    model_id: Either the integer model ID to parse, or one of ModelID.FIRST to
      parse the first model, or ModelID.ALL to parse all models.

  Returns:
    A Structure representation of the mmCIF object.
  """
  str_model_id = _get_str_model_id(cif=mmcif_object, model_id=model_id)
  header = _get_mmcif_header(
      mmcif_object, fix_mse=fix_mse_residues, fix_unknown_dna=fix_unknown_dna
  )

  chains, residues, atoms = get_tables(
      cif=mmcif_object,
      fix_mse_residues=fix_mse_residues,
      fix_arginines=fix_arginines,
      fix_unknown_dna=fix_unknown_dna,
      include_water=include_water,
      include_other=include_other,
      model_id=str_model_id,
  )

  if include_bonds and atoms.size > 0:
    # NB: parsing the atom table before the bonds table allows for a more
    # informative error message when dealing with bad multi-model mmCIFs.
    # Also always use a specific model ID, even when parsing all models.
    if str_model_id == '':  # pylint: disable=g-explicit-bool-comparison
      bonds_model_id = _get_first_model_id(mmcif_object)
    else:
      bonds_model_id = str_model_id

    bonds_table = _parse_bonds(
        mmcif_object,
        atom_key=atoms.key,
        model_id=bonds_model_id,
    )
  else:
    bonds_table = bonds.Bonds.make_empty()

  return structure.Structure(
      name=name if name is not None else header.name,
      resolution=header.resolution,
      release_date=header.release_date,
      structure_method=header.structure_method,
      bioassembly_data=header.bioassembly_data,
      chemical_components_data=header.chemical_components_data,
      bonds=bonds_table,
      chains=chains,
      residues=residues,
      atoms=atoms,
  )


def from_mmcif(
    mmcif_string: str | bytes,
    *,
    name: str | None = None,
    fix_mse_residues: bool = False,
    fix_arginines: bool = False,
    fix_unknown_dna: bool = False,
    include_water: bool = False,
    include_other: bool = False,
    include_bonds: bool = False,
    model_id: int | ModelID = ModelID.FIRST,
) -> structure.Structure:
  """Construct a Structure from a mmCIF string.

  Args:
    mmcif_string: The string contents of an mmCIF file.
    name: Optional name for the structure. If not provided, the name will be
      taken from the mmCIF data_ field.
    fix_mse_residues: If True, selenium atom sites (SE) in selenomethionine
      (MSE) residues will be changed to sulphur atom sites (SD). This is because
      methionine (MET) residues are often replaced with MSE to aid X-Ray
      crystallography. If False, the SE MSE atom sites won't be modified.
    fix_arginines: If True, NH1 and NH2 in arginine will be swapped if needed so
      that NH1 is always closer to CD than NH2. If False, no atom sites in
      arginine will be touched. Note that HH11, HH12, HH21, HH22 are fixed too.
    fix_unknown_dna: If True, residues with name N in DNA chains will have their
      res_name replaced with DN. Atoms are not changed.
    include_water: If True, water (HOH) molecules will be parsed. Water
      molecules may be grouped into chains, where number of residues > 1. Water
      molecules are usually grouped into chains but do not necessarily all share
      the same chain ID.
    include_other: If True, all other atoms that are not included by any of the
      above parameters will be included. This covers e.g. "polypeptide(D)" and
      "macrolide" entities, as well as all other non-standard types.
    include_bonds: If True, bond information will be parsed from the mmCIF and
      stored in the Structure.
    model_id: Either the integer model ID to parse, or one of ModelID.FIRST to
      parse the first model, or ModelID.ALL to parse all models.

  Returns:
    A Structure representation of the mmCIF string.
  """
  mmcif_object = mmcif.from_string(mmcif_string)

  return from_parsed_mmcif(
      mmcif_object,
      name=name,
      fix_mse_residues=fix_mse_residues,
      fix_arginines=fix_arginines,
      fix_unknown_dna=fix_unknown_dna,
      include_water=include_water,
      include_other=include_other,
      include_bonds=include_bonds,
      model_id=model_id,
  )


def from_res_arrays(atom_mask: np.ndarray, **kwargs) -> structure.Structure:
  """Returns Structure created from from arrays with a residue dimension.

  All unset fields are filled with defaults (e.g. 1.0 for occupancy) or
  unset/unknown values (e.g. UNK for residue type, or '.' for atom element).

  Args:
    atom_mask: A array with shape (num_res, num_atom). This is used to decide
      which atoms in the atom dimension are present in a given residue. Present
      atoms should have a nonzero value, e.g. 1.0 or True.
    **kwargs: A mapping from field name to values. For all array-valued fields
      these arrays must have a dimension of length num_res. Chain and residue
      fields should have this as their only dimension and atom fields should be
      shaped (num_res, num_atom). Coordinate fields may also have arbitrary
      leading dimensions (they must be the same across all coordinate fields).
      See structure.{CHAIN,RESIDUE,ATOM}_FIELDS for a list of allowed fields.
  """
  num_res, num_atom = atom_mask.shape
  included_indices = np.flatnonzero(atom_mask)

  array_fields = (
      structure.CHAIN_FIELDS.keys()
      | structure.RESIDUE_FIELDS.keys()
      | structure.ATOM_FIELDS.keys()
  )
  initializer_kwargs = {}
  fields = {}
  for k, val in kwargs.items():
    if k not in array_fields:
      # The kwarg key isn't an array field name. Such kwargs are forwarded as-is
      # to the constructor. They are expected to be global fields (e.g. name).
      # Other values will raise an error when the constructor is called.
      if k in structure.TABLE_FIELDS:
        raise ValueError(f'Table fields must not be set. Got {k}.')
      initializer_kwargs[k] = val
      continue
    elif val is None:
      raise ValueError(f'{k} must be non-None.')

    if not isinstance(val, np.ndarray):
      raise TypeError(f'Value for {k} must be a NumPy array. Got {type(val)}.')
    if k in structure.CHAIN_FIELDS or k in structure.RESIDUE_FIELDS:
      if val.shape != (num_res,):
        raise ValueError(
            f'{k} must have shape ({num_res=},). Got {val.shape=}.'
        )
      # Do not reshape the chain/residue arrays, they have the shape we need.
      fields[k] = val
    else:
      assert k in structure.ATOM_FIELDS
      if val.shape[-2:] != (num_res, num_atom):
        raise ValueError(
            f'{k} must have final two dimensions of length '
            f'{(num_res, num_atom)=}. Got {val.shape=}.'
        )
      leading_dims = val.shape[:-2]
      flat_val = val.reshape(leading_dims + (-1,), order='C')
      masked_val = flat_val[..., included_indices]
      fields[k] = masked_val

  # Get chain IDs or assume this is a single-chain structure.
  chain_id = kwargs.get('chain_id', np.array(['A'] * num_res, dtype=object))
  # Find chain starts in res-sized arrays, use these to make chain-sized arrays.
  chain_start = np.concatenate(
      ([0], np.where(chain_id[1:] != chain_id[:-1])[0] + 1)
  )
  if len(set(chain_id)) != len(chain_start):
    raise ValueError(f'Chain IDs must be contiguous, but got {chain_id}')

  chain_lengths = np.diff(chain_start, append=len(chain_id))
  chain_key = np.repeat(np.arange(len(chain_start)), chain_lengths)

  chain_entity_id = fields.get('chain_entity_id')
  if chain_entity_id is not None:
    entity_id = chain_entity_id[chain_start]
  else:
    entity_id = np.array(
        [str(mmcif.str_id_to_int_id(cid)) for cid in chain_id[chain_start]],
        dtype=object,
    )
  chain_str_empty = np.full((num_res,), '.', dtype=object)
  chains_table = structure_tables.Chains(
      key=chain_key[chain_start],
      id=chain_id[chain_start],
      type=fields.get('chain_type', chain_str_empty)[chain_start],
      auth_asym_id=fields.get('chain_auth_asym_id', chain_id)[chain_start],
      entity_id=entity_id,
      entity_desc=fields.get('chain_entity_desc', chain_str_empty)[chain_start],
  )

  # Since all arrays are residue-shaped, we can use them directly.
  res_key = np.arange(num_res, dtype=np.int64)
  res_id = fields.get('res_id', res_key + 1).astype(np.int32)
  residues_table = structure_tables.Residues(
      key=res_key,
      chain_key=chain_key,
      id=res_id,
      name=fields.get('res_name', np.full(num_res, 'UNK', dtype=object)),
      auth_seq_id=fields.get(
          'res_auth_seq_id', np.char.mod('%d', res_id).astype(object)
      ),
      insertion_code=fields.get(
          'res_insertion_code', np.full(num_res, '?', dtype=object)
      ),
  )

  # The atom-sized arrays have already been masked and reshaped.
  num_atoms_per_res = np.sum(atom_mask, axis=1, dtype=np.int32)
  num_atoms_total = np.sum(num_atoms_per_res, dtype=np.int32)
  # Structure is immutable, so use the same array multiple times to save RAM.
  atom_str_empty = np.full(num_atoms_total, '.', dtype=object)
  atom_float32_zeros = np.zeros(num_atoms_total, dtype=np.float32)
  atom_float32_ones = np.ones(num_atoms_total, dtype=np.float32)
  atoms_table = structure_tables.Atoms(
      key=np.arange(num_atoms_total, dtype=np.int64),
      chain_key=np.repeat(chain_key, num_atoms_per_res),
      res_key=np.repeat(res_key, num_atoms_per_res),
      name=fields.get('atom_name', atom_str_empty),
      element=fields.get('atom_element', atom_str_empty),
      x=fields.get('atom_x', atom_float32_zeros),
      y=fields.get('atom_y', atom_float32_zeros),
      z=fields.get('atom_z', atom_float32_zeros),
      b_factor=fields.get('atom_b_factor', atom_float32_zeros),
      occupancy=fields.get('atom_occupancy', atom_float32_ones),
  )

  return structure.Structure(
      chains=chains_table,
      residues=residues_table,
      atoms=atoms_table,
      bonds=structure_tables.Bonds.make_empty(),  # Currently not set.
      **initializer_kwargs,
  )


def expand_sequence(
    sequence: str, chain_type: str, sequence_format: SequenceFormat
) -> Sequence[str]:
  """Returns full residue names based on a sequence string.

  Args:
    sequence: A string representing the sequence.
    chain_type: The chain type of the sequence.
    sequence_format: The format of the sequence argument.
  """
  match sequence_format:
    case SequenceFormat.FASTA:
      if not all(c.isalpha() for c in sequence):
        raise ValueError(f'Sequence "{sequence}" has non-alphabetic characters')
      match chain_type:
        case mmcif_names.PROTEIN_CHAIN:
          res_name_map = residue_names.PROTEIN_COMMON_ONE_TO_THREE
          default_res_name = residue_names.UNK
        case mmcif_names.RNA_CHAIN:
          res_name_map = {r: r for r in residue_names.RNA_TYPES}
          default_res_name = residue_names.UNK_RNA
        case mmcif_names.DNA_CHAIN:
          res_name_map = residue_names.DNA_COMMON_ONE_TO_TWO
          default_res_name = residue_names.UNK_DNA
        case _:
          raise ValueError(f'{chain_type=} not supported for FASTA format.')
      return [
          res_name_map.get(one_letter_res, default_res_name)
          for one_letter_res in sequence
      ]
    case SequenceFormat.CCD_CODES:
      return sequence.strip('()').split(')(')
    case SequenceFormat.LIGAND_SMILES:
      ligand_id, _ = sequence.split(':', maxsplit=1)
      return [ligand_id]


def from_sequences_and_bonds(
    *,
    sequences: Sequence[str],
    chain_types: Sequence[str],
    sequence_formats: Sequence[SequenceFormat],
    bonded_atom_pairs: Sequence[tuple[BondAtomId, BondAtomId]] | None,
    ccd: chemical_components.Ccd,
    chain_ids: Sequence[str] | None = None,
    name: str = 'from_sequences_and_bonds',
    bond_type: str | None = None,
    **constructor_args,
) -> structure.Structure:
  """Returns a minimal structure for the input sequences and bonds.

  The returned structure will have at least one atom per residue. If the
  residue has any bonded atoms, according to `bonded_atom_pairs`, then
  all (and only) those atoms will be present for that residue. If the residue
  is not involved in any bond then an arbitrary atom will be created.

  Args:
    sequences: A sequence of strings, each one representing a single chain.
    chain_types: The types of each chain, e.g. polypeptide(L). The n-th element
      describes the n-th sequence in `sequences`.
    sequence_formats: The format of each sequence. The n-th element describes
      the n-th sequence in `sequences`.
    bonded_atom_pairs: A sequence of bonded atom pairs. Each atom is described
      as a tuple of (chain_index, res_index, atom_name), where the first two
      values are 0-based indices. The chain_index is the index of the chain in
      the `sequences` argument, and the res_index is the index of the residue in
      that sequence. The atom_name is the name of the atom in the residue, e.g.
      CA. If the atom is not found in the standard atoms for that residue
      (according to the CCD) then an error is raised.
    ccd: The chemical components dictionary.
    chain_ids: A sequence of chain IDs, one for each chain in `sequences`. If
      not provided, then the chain IDs will be generated automatically based on
      sequence indices.
    name: A name for the returned structure.
    bond_type: This type will be used for all bonds in the structure, where type
      follows PDB scheme, e.g. unknown (?), hydrog, metalc, covale, disulf.
    **constructor_args: These arguments are passed directly to the
      structure.Structure constructor.
  """
  chain_id = []
  chain_type = []
  chain_res_count = []
  res_id = []
  res_name = []
  res_atom_count = []
  atom_name = []
  atom_element = []
  chem_comp = {}

  num_bonds = len(bonded_atom_pairs or ())
  from_atom_key = np.full((num_bonds,), -1, dtype=np.int64)
  dest_atom_key = np.full((num_bonds,), -1, dtype=np.int64)

  # Create map (chain_i, res_i) -> {atom_name -> (from_idxs dest_idxs)}.
  # This allows quick lookup of whether a residue has any bonded atoms, and
  # which bonds those atoms participate in.
  bond_lookup = _create_bond_lookup(bonded_atom_pairs or ())

  current_atom_key = 0
  for chain_i, (sequence, curr_chain_type, sequence_format) in enumerate(
      zip(sequences, chain_types, sequence_formats, strict=True)
  ):
    if chain_ids is not None:
      current_chain_id = chain_ids[chain_i]
    else:
      current_chain_id = mmcif.int_id_to_str_id(chain_i + 1)
    num_chain_residues = 0
    for res_i, full_res_name in enumerate(
        expand_sequence(sequence, curr_chain_type, sequence_format)
    ):
      current_res_id = res_i + 1
      num_res_atoms = 0

      # Look for bonded atoms in the bond lookup and if any are found, add
      # their atom keys to the bond atom_key columns.
      if bond_indices_by_atom_name := bond_lookup.get((chain_i, res_i)):
        comp_atoms = None
        if sequence_format != SequenceFormat.LIGAND_SMILES:
          comp_atoms = set(ccd.get(full_res_name)['_chem_comp_atom.atom_id'])
        for bond_atom_name, bond_indices in bond_indices_by_atom_name.items():
          if comp_atoms is not None and bond_atom_name not in comp_atoms:
            raise ValueError(
                f'Bonded atom "{bond_atom_name}" was not found in the list of'
                f' atoms of the chemical component {full_res_name}. Valid atom'
                f' names for {full_res_name} are: {sorted(comp_atoms)}.'
                ' This is likely caused by an invalid atom name in the bonded'
                f' atom (chain_id={current_chain_id}, res_id={current_res_id},'
                f' atom_name={bond_atom_name}) specified in `bondedAtomPairs`'
                ' in the input JSON.'
            )
          atom_name.append(bond_atom_name)
          atom_element.append(
              _get_atom_element(
                  ccd=ccd, res_name=full_res_name, atom_name=bond_atom_name
              )
          )
          for from_bond_i in bond_indices.from_indices:
            from_atom_key[from_bond_i] = current_atom_key
          for dest_bond_i in bond_indices.dest_indices:
            dest_atom_key[dest_bond_i] = current_atom_key
          current_atom_key += 1
          num_res_atoms += 1
      else:
        # If this residue has no bonded atoms then we need to add one atom
        # like in from_sequences.
        assert num_res_atoms == 0
        rep_atom_name, rep_atom_element = _get_representative_atom(
            ccd=ccd,
            res_name=full_res_name,
            chain_type=curr_chain_type,
            sequence_format=sequence_format,
        )
        atom_name.append(rep_atom_name)
        atom_element.append(rep_atom_element)
        num_res_atoms += 1
        current_atom_key += 1

      if sequence_format == SequenceFormat.LIGAND_SMILES:
        # Sequence expect to be in the format <ligand_id>:<ligand_smiles>,
        # which always corresponds to a single-residue chain.
        ligand_id, ligand_smiles = sequence.split(':', maxsplit=1)
        if ccd.get(ligand_id) is not None:
          raise ValueError(
              f'Ligand name {ligand_id} is in CCD - it is not supported to give'
              ' ligands created from SMILES the same name as CCD components.'
          )
        # We need to provide additional chemical components metadata for
        # ligands specified via SMILES strings since they might not be in CCD.
        _add_ligand_to_chem_comp(chem_comp, ligand_id, ligand_smiles)

      assert num_res_atoms >= 1
      res_atom_count.append(num_res_atoms)
      num_chain_residues += 1
      res_id.append(current_res_id)
      res_name.append(full_res_name)

    chain_id.append(current_chain_id)
    chain_type.append(curr_chain_type)
    chain_res_count.append(num_chain_residues)

  chem_comp_data = struc_chem_comps.ChemicalComponentsData(chem_comp)
  chem_comp_data = struc_chem_comps.populate_missing_ccd_data(
      ccd=ccd,
      chemical_components_data=chem_comp_data,
      chemical_component_ids=set(res_name),
  )

  if bonded_atom_pairs is not None:
    unknown_bond_col = np.full((num_bonds,), '?', dtype=object)
    if bond_type is None:
      bond_type_col = unknown_bond_col
    else:
      bond_type_col = np.full((num_bonds,), bond_type, dtype=object)
    bonds_table = bonds.Bonds(
        key=np.arange(num_bonds, dtype=np.int64),
        type=bond_type_col,
        role=unknown_bond_col,
        from_atom_key=from_atom_key,
        dest_atom_key=dest_atom_key,
    )
  else:
    bonds_table = structure_tables.Bonds.make_empty()

  chain_key = np.arange(len(sequences), dtype=np.int64)  # 1 chain per sequence.
  chain_id = np.array(chain_id, dtype=object)
  chains_table = structure_tables.Chains(
      key=chain_key,
      id=chain_id,
      type=np.array(chain_type, dtype=object),
      auth_asym_id=chain_id,
      entity_id=np.char.mod('%d', chain_key + 1).astype(object),
      entity_desc=np.array(['.'] * len(chain_key), dtype=object),
  )

  res_key = np.arange(len(res_name), dtype=np.int64)
  res_chain_key = np.repeat(chain_key, chain_res_count)
  residues_table = structure_tables.Residues(
      key=res_key,
      chain_key=res_chain_key,
      id=np.array(res_id, dtype=np.int32),
      name=np.array(res_name, dtype=object),
      auth_seq_id=np.char.mod('%d', res_id).astype(object),
      insertion_code=np.full(len(res_name), '?', dtype=object),
  )

  num_atoms = current_atom_key
  atom_float32_zeros = np.zeros(num_atoms, dtype=np.float32)
  atoms_table = structure_tables.Atoms(
      key=np.arange(num_atoms, dtype=np.int64),
      chain_key=np.repeat(res_chain_key, res_atom_count),
      res_key=np.repeat(res_key, res_atom_count),
      name=np.array(atom_name, dtype=object),
      element=np.array(atom_element, dtype=object),
      x=atom_float32_zeros,
      y=atom_float32_zeros,
      z=atom_float32_zeros,
      b_factor=atom_float32_zeros,
      occupancy=np.ones(num_atoms, np.float32),
  )

  return structure.Structure(
      name=name,
      atoms=atoms_table,
      residues=residues_table,
      chains=chains_table,
      bonds=bonds_table,
      chemical_components_data=chem_comp_data,
      **constructor_args,
  )


class _ChainResBuilder:
  """Class for incrementally building chain and residue tables."""

  def __init__(
      self,
      *,
      chain_key_by_chain_id: Mapping[str, int],
      entity_id_by_chain_id: Mapping[str, str],
      chain_type_by_entity_id: Mapping[str, str],
      entity_desc_by_entity_id: Mapping[str, str],
      fix_mse_residues: bool,
      fix_unknown_dna: bool,
  ):
    # Len: num_chains.
    self.chain_key = []
    self.chain_id = []
    self.chain_type = []
    self.chain_auth_asym_id = []
    self.chain_entity_id = []
    self.chain_entity_desc = []

    # Len: num_residues.
    self.res_key = []
    self.res_chain_key = []
    self.res_id = []
    self.res_name = []
    self.res_auth_seq_id = []
    self.res_insertion_code = []

    self.chain_key_by_chain_id = chain_key_by_chain_id
    self.entity_id_by_chain_id = entity_id_by_chain_id
    self.chain_type_by_entity_id = chain_type_by_entity_id
    self.entity_desc_by_entity_id = entity_desc_by_entity_id
    self.key_for_res: dict[tuple[str, str, str, str], int] = {}

    self._fix_mse_residues = fix_mse_residues
    self._fix_unknown_dna = fix_unknown_dna

  def add_residues(
      self,
      *,
      chain_ids: np.ndarray,
      chain_auth_asym_ids: np.ndarray,
      res_ids: np.ndarray,
      res_names: np.ndarray,
      res_auth_seq_ids: np.ndarray,
      res_ins_codes: np.ndarray,
  ):
    """Adds a residue (and its chain) to the tables."""
    # Create chain table data.
    if chain_ids.size == 0:
      return

    chain_ids_with_prev = np.concatenate(
        (([self.chain_id[-1] if self.chain_id else None], chain_ids))
    )
    chain_change_mask = chain_ids_with_prev[:-1] != chain_ids_with_prev[1:]
    chain_change_ids = chain_ids[chain_change_mask]
    chain_keys = string_array.remap(
        chain_change_ids, self.chain_key_by_chain_id, inplace=False
    )
    self.chain_key.extend(chain_keys)
    self.chain_id.extend(chain_change_ids)
    self.chain_auth_asym_id.extend(chain_auth_asym_ids[chain_change_mask])
    chain_entity_id = string_array.remap(
        chain_change_ids, self.entity_id_by_chain_id, inplace=False
    )
    self.chain_entity_id.extend(chain_entity_id)
    chain_type = string_array.remap(
        chain_entity_id, self.chain_type_by_entity_id, inplace=False
    )
    self.chain_type.extend(chain_type)
    chain_entity_desc = string_array.remap(
        chain_entity_id, self.entity_desc_by_entity_id, inplace=False
    )
    self.chain_entity_desc.extend(chain_entity_desc)

    # Create residue table data.
    num_prev_res = len(self.res_id)
    res_keys = np.arange(num_prev_res, num_prev_res + len(res_ids))
    res_iter = zip(
        chain_ids,
        res_auth_seq_ids,
        res_names,
        res_ins_codes,
        strict=True,
    )
    key_for_res_update = {
        res_unique_id: res_key
        for res_key, res_unique_id in enumerate(res_iter, num_prev_res)
    }
    self.key_for_res.update(key_for_res_update)
    self.res_key.extend(res_keys)
    self.res_chain_key.extend(
        string_array.remap(chain_ids, self.chain_key_by_chain_id, inplace=False)
    )
    self.res_id.extend(res_ids)
    self.res_name.extend(res_names)
    self.res_auth_seq_id.extend(res_auth_seq_ids)
    self.res_insertion_code.extend(res_ins_codes)

  def make_chains_table(self) -> structure_tables.Chains:
    """Returns the Structure chains table."""
    chain_key = np.array(self.chain_key, dtype=np.int64)
    if not np.all(chain_key[:-1] <= chain_key[1:]):
      # If the order is inconsistent with the atoms table, sort so that it is.
      order = np.argsort(self.chain_key, kind='stable')
      return structure_tables.Chains(
          key=chain_key[order],
          id=np.array(self.chain_id, dtype=object)[order],
          type=np.array(self.chain_type, dtype=object)[order],
          auth_asym_id=np.array(self.chain_auth_asym_id, dtype=object)[order],
          entity_id=np.array(self.chain_entity_id, dtype=object)[order],
          entity_desc=np.array(self.chain_entity_desc, dtype=object)[order],
      )
    return structure_tables.Chains(
        key=chain_key,
        id=np.array(self.chain_id, dtype=object),
        type=np.array(self.chain_type, dtype=object),
        auth_asym_id=np.array(self.chain_auth_asym_id, dtype=object),
        entity_id=np.array(self.chain_entity_id, dtype=object),
        entity_desc=np.array(self.chain_entity_desc, dtype=object),
    )

  def make_residues_table(self) -> structure_tables.Residues:
    """Returns the Structure residues table."""
    res_name = np.array(self.res_name, dtype=object)
    res_chain_key = np.array(self.res_chain_key, dtype=np.int64)

    if self._fix_mse_residues:
      string_array.remap(res_name, mapping={'MSE': 'MET'}, inplace=True)

    if self._fix_unknown_dna:
      # Remap residues from N -> DN in DNA chains only.
      dna_chain_mask = (
          np.array(self.chain_type, dtype=object) == mmcif_names.DNA_CHAIN
      )
      dna_chain_key = np.array(self.chain_key, dtype=object)[dna_chain_mask]
      res_name[(res_name == 'N') & np.isin(res_chain_key, dna_chain_key)] = 'DN'

    if not np.all(res_chain_key[:-1] <= res_chain_key[1:]):
      # If the order is inconsistent with the atoms table, sort so that it is.
      order = np.argsort(res_chain_key, kind='stable')
      return structure_tables.Residues(
          key=np.array(self.res_key, dtype=np.int64)[order],
          chain_key=res_chain_key[order],
          id=np.array(self.res_id, dtype=np.int32)[order],
          name=res_name[order],
          auth_seq_id=np.array(self.res_auth_seq_id, dtype=object)[order],
          insertion_code=np.array(self.res_insertion_code, dtype=object)[order],
      )
    return structure_tables.Residues(
        key=np.array(self.res_key, dtype=np.int64),
        chain_key=res_chain_key,
        id=np.array(self.res_id, dtype=np.int32),
        name=res_name,
        auth_seq_id=np.array(self.res_auth_seq_id, dtype=object),
        insertion_code=np.array(self.res_insertion_code, dtype=object),
    )


def _get_string_array_default(cif: mmcif.Mmcif, key: str, default: list[str]):
  try:
    return cif.get_array(key, dtype=object)
  except KeyError:
    return default


def _generate_required_tables_if_missing(
    cif: mmcif.Mmcif,
) -> Mapping[str, Sequence[str]]:
  """Generates all required tables and columns if missing."""
  update = {}

  atom_site_entities = _get_string_array_default(
      cif, '_atom_site.label_entity_id', []
  )

  # OpenMM produces files that don't have any of the tables and also have
  # _atom_site.label_entity_id set to '?' for all atoms. We infer the entities
  # based on the _atom_site.label_asym_id column. We start with cheaper O(1)
  # checks to prevent running the expensive O(n) check on most files.
  if (
      len(atom_site_entities) > 0  # pylint: disable=g-explicit-length-test
      and '_entity.id' not in cif  # Ignore if the _entity table exists.
      and atom_site_entities[0] == '?'  # Cheap check.
      and set(atom_site_entities) == {'?'}  # Expensive check.
  ):
    label_asym_ids = cif.get_array('_atom_site.label_asym_id', dtype=object)
    atom_site_entities = [
        str(mmcif.str_id_to_int_id(cid)) for cid in label_asym_ids
    ]
    # Update _atom_site.label_entity_id to be consistent with the new tables.
    update['_atom_site.label_entity_id'] = atom_site_entities

  # Check table existence by checking the presence of its primary key.
  if '_struct_asym.id' not in cif:
    # Infer the _struct_asym table using the _atom_site table.
    asym_ids = _get_string_array_default(cif, '_atom_site.label_asym_id', [])

    if len(atom_site_entities) == 0 or len(asym_ids) == 0:  # pylint: disable=g-explicit-length-test
      raise ValueError(
          'Could not parse an mmCIF with no _struct_asym table and also no '
          '_atom_site.label_entity_id or _atom_site.label_asym_id columns.'
      )

    # Deduplicate, but keep the order intact - dict.fromkeys maintains order.
    entity_id_chain_id_pairs = list(
        dict.fromkeys(zip(atom_site_entities, asym_ids, strict=True))
    )
    update['_struct_asym.entity_id'] = [e for e, _ in entity_id_chain_id_pairs]
    update['_struct_asym.id'] = [c for _, c in entity_id_chain_id_pairs]

  if '_entity.id' not in cif:
    # Infer the _entity_poly and _entity tables using the _atom_site table.
    residues = _get_string_array_default(cif, '_atom_site.label_comp_id', [])
    group_pdb = _get_string_array_default(cif, '_atom_site.group_PDB', [])
    if '_atom_site.label_entity_id' in cif:
      entities = atom_site_entities
    else:
      # If _atom_site.label_entity_id not set, use the asym_id -> entity_id map.
      asym_to_entity = dict(
          zip(
              cif['_struct_asym.id'], cif['_struct_asym.entity_id'], strict=True
          )
      )
      entities = string_array.remap(
          cif.get_array('_atom_site.label_asym_id', dtype=object),
          mapping=asym_to_entity,
      )

    entity_ids = []
    entity_types = []
    entity_poly_entity_ids = []
    entity_poly_types = []
    entity_poly_table_missing = '_entity_poly.entity_id' not in cif
    for entity_id, group in itertools.groupby(
        zip(entities, residues, group_pdb, strict=True), key=lambda e: e[0]
    ):
      _, entity_residues, entity_group_pdb = zip(*group, strict=True)
      entity_type = _guess_entity_type(
          chain_residues=entity_residues, atom_types=entity_group_pdb
      )
      entity_ids.append(entity_id)
      entity_types.append(entity_type)

      if entity_poly_table_missing and entity_type == mmcif_names.POLYMER_CHAIN:
        polymer_type = mmcif_names.guess_polymer_type(entity_residues)
        entity_poly_entity_ids.append(entity_id)
        entity_poly_types.append(polymer_type)

    update['_entity.id'] = entity_ids
    update['_entity.type'] = entity_types
    if entity_poly_table_missing:
      update['_entity_poly.entity_id'] = entity_poly_entity_ids
      update['_entity_poly.type'] = entity_poly_types

  if '_atom_site.type_symbol' not in cif:
    update['_atom_site.type_symbol'] = mmcif.get_or_infer_type_symbol(cif)

  return update


def _maybe_add_missing_scheme_tables(
    cif: mmcif.Mmcif,
    res_starts: Sequence[int],
    label_asym_ids: np.ndarray,
    label_seq_ids: np.ndarray,
    label_comp_ids: np.ndarray,
    auth_seq_ids: np.ndarray,
    pdb_ins_codes: np.ndarray,
) -> Mapping[str, Sequence[str]]:
  """If missing, infers the scheme tables from the _atom_site table."""
  update = {}

  required_poly_seq_scheme_cols = (
      '_pdbx_poly_seq_scheme.asym_id',
      '_pdbx_poly_seq_scheme.pdb_seq_num',
      '_pdbx_poly_seq_scheme.pdb_ins_code',
      '_pdbx_poly_seq_scheme.seq_id',
      '_pdbx_poly_seq_scheme.mon_id',
      '_pdbx_poly_seq_scheme.pdb_strand_id',
  )
  if not all(col in cif for col in required_poly_seq_scheme_cols):
    # Create a mask for atoms where each polymer residue start.
    entity_id_by_chain_id = dict(
        zip(cif['_struct_asym.id'], cif['_struct_asym.entity_id'], strict=True)
    )
    chain_type_by_entity_id = dict(
        zip(cif['_entity.id'], cif['_entity.type'], strict=True)
    )
    # Remap asym ID -> entity ID.
    label_entity_id = string_array.remap(
        label_asym_ids, mapping=entity_id_by_chain_id, inplace=False
    )
    # Remap entity ID -> chain type.
    chain_type = string_array.remap(
        label_entity_id, mapping=chain_type_by_entity_id, inplace=False
    )
    res_mask = np.zeros_like(label_seq_ids, dtype=bool)
    res_mask[res_starts] = True
    res_mask &= chain_type == mmcif_names.POLYMER_CHAIN

    entity_poly_seq_cols = (
        '_entity_poly_seq.entity_id',
        '_entity_poly_seq.num',
        '_entity_poly_seq.mon_id',
    )
    if all(col in cif for col in entity_poly_seq_cols):
      # Use _entity_poly_seq if available.
      poly_seq_num = cif.get_array('_entity_poly_seq.num', dtype=object)
      poly_seq_mon_id = cif.get_array('_entity_poly_seq.mon_id', dtype=object)
      poly_seq_entity_id = cif.get_array(
          '_entity_poly_seq.entity_id', dtype=object
      )
      # We have to add the entity ID to the residue ID because multiple residues
      # can share the same ID. This also allows using string_array.remap.
      label_seq_id_to_auth_seq_id = dict(
          zip(
              np.char.add(label_entity_id[res_mask], label_seq_ids[res_mask]),
              auth_seq_ids[res_mask],
              strict=True,
          )
      )
      scheme_pdb_seq_num = string_array.remap(
          np.char.add(poly_seq_entity_id, poly_seq_num),
          mapping=label_seq_id_to_auth_seq_id,
          default_value='.',
      )
      label_seq_id_to_ins_code = dict(
          zip(label_seq_ids[res_mask], pdb_ins_codes[res_mask], strict=True)
      )
      scheme_pdb_ins_code = string_array.remap(
          poly_seq_num, mapping=label_seq_id_to_ins_code, default_value='.'
      )

      # The _entity_poly_seq table is entity-based, while _pdbx_poly_seq_scheme
      # is chain-based. A single entity could mean multiple chains (asym_ids),
      # we therefore need to replicate each entity for all of the chains.
      scheme_asym_id = []
      select = []
      indices = np.arange(len(poly_seq_entity_id), dtype=np.int32)
      for asym_id, entity_id in zip(
          cif['_struct_asym.id'], cif['_struct_asym.entity_id'], strict=True
      ):
        entity_mask = poly_seq_entity_id == entity_id
        select.extend(indices[entity_mask])
        scheme_asym_id.extend([asym_id] * sum(entity_mask))

      scheme_pdb_strand_id = string_array.remap(
          np.array(scheme_asym_id, dtype=object),
          mapping=mmcif.get_internal_to_author_chain_id_map(cif),
          inplace=False,
      )

      update['_pdbx_poly_seq_scheme.asym_id'] = scheme_asym_id
      update['_pdbx_poly_seq_scheme.pdb_strand_id'] = scheme_pdb_strand_id
      update['_pdbx_poly_seq_scheme.pdb_seq_num'] = scheme_pdb_seq_num[select]
      update['_pdbx_poly_seq_scheme.pdb_ins_code'] = scheme_pdb_ins_code[select]
      update['_pdbx_poly_seq_scheme.seq_id'] = poly_seq_num[select]
      update['_pdbx_poly_seq_scheme.mon_id'] = poly_seq_mon_id[select]
    else:
      # _entity_poly_seq not available, fallback to _atom_site.
      res_asym_ids = label_asym_ids[res_mask]
      res_strand_ids = string_array.remap(
          array=res_asym_ids,
          mapping=mmcif.get_internal_to_author_chain_id_map(cif),
          inplace=False,
      )
      update['_pdbx_poly_seq_scheme.asym_id'] = res_asym_ids
      update['_pdbx_poly_seq_scheme.pdb_strand_id'] = res_strand_ids
      update['_pdbx_poly_seq_scheme.pdb_seq_num'] = auth_seq_ids[res_mask]
      update['_pdbx_poly_seq_scheme.pdb_ins_code'] = pdb_ins_codes[res_mask]
      update['_pdbx_poly_seq_scheme.seq_id'] = label_seq_ids[res_mask]
      update['_pdbx_poly_seq_scheme.mon_id'] = label_comp_ids[res_mask]

  required_nonpoly_scheme_cols = (
      '_pdbx_nonpoly_scheme.mon_id',
      '_pdbx_nonpoly_scheme.asym_id',
      '_pdbx_nonpoly_scheme.pdb_seq_num',
      '_pdbx_nonpoly_scheme.pdb_ins_code',
  )
  required_branch_scheme_cols = (
      '_pdbx_branch_scheme.mon_id',
      '_pdbx_branch_scheme.asym_id',
      '_pdbx_branch_scheme.pdb_seq_num',
  )

  # Generate _pdbx_nonpoly_scheme only if both tables are missing.
  if not (
      all(col in cif for col in required_nonpoly_scheme_cols)
      or all(col in cif for col in required_branch_scheme_cols)
  ):
    # To be strictly semantically correct, multi-residue ligands should be
    # written in _pdbx_branch_scheme. However, Structure parsing handles
    # correctly multi-residue ligands in _pdbx_nonpoly_scheme and the tables
    # constructed here live only while parsing, hence this is unnecessary.
    entity_id_by_chain_id = dict(
        zip(cif['_struct_asym.id'], cif['_struct_asym.entity_id'], strict=True)
    )
    chain_type_by_entity_id = dict(
        zip(cif['_entity.id'], cif['_entity.type'], strict=True)
    )
    # Remap asym ID -> entity ID.
    chain_type = string_array.remap(
        label_asym_ids, mapping=entity_id_by_chain_id, inplace=False
    )
    # Remap entity ID -> chain type.
    string_array.remap(
        chain_type, mapping=chain_type_by_entity_id, inplace=True
    )
    res_mask = np.zeros_like(label_seq_ids, dtype=bool)
    res_mask[res_starts] = True
    res_mask &= chain_type != mmcif_names.POLYMER_CHAIN

    if not np.any(res_mask):
      return update  # Shortcut: no non-polymer residues.

    ins_codes = string_array.remap(
        pdb_ins_codes[res_mask], mapping={'?': '.'}, inplace=False
    )

    update['_pdbx_nonpoly_scheme.asym_id'] = label_asym_ids[res_mask]
    update['_pdbx_nonpoly_scheme.pdb_seq_num'] = auth_seq_ids[res_mask]
    update['_pdbx_nonpoly_scheme.pdb_ins_code'] = ins_codes
    update['_pdbx_nonpoly_scheme.mon_id'] = label_comp_ids[res_mask]

  return update


def _get_chain_key_by_chain_id(
    resolved_chain_ids: np.ndarray, struct_asym_chain_ids: np.ndarray
) -> Mapping[str, int]:
  """Returns chain key for each chain ID respecting resolved chain ordering."""
  # Check that all chain IDs found in the (potentially filtered) _atom_site
  # table are present in the _struct_asym table.
  unique_resolved_chain_ids = set(resolved_chain_ids)
  if not unique_resolved_chain_ids.issubset(set(struct_asym_chain_ids)):
    unique_resolved_chain_ids = sorted(unique_resolved_chain_ids)
    unique_struct_asym_chain_ids = sorted(set(struct_asym_chain_ids))
    raise ValueError(
        'Bad mmCIF: chain IDs in _atom_site.label_asym_id '
        f'{unique_resolved_chain_ids} is not a subset of chain IDs in '
        f'_struct_asym.id {unique_struct_asym_chain_ids}.'
    )

  resolved_mask = string_array.isin(
      struct_asym_chain_ids, unique_resolved_chain_ids
  )
  # For all resolved chains, use the _atom_site order they appear in. E.g.
  # resolved_chain_ids     = [B A   E D F]
  # struct_asym_chain_ids  = [A B C D E F]
  # consistent_chain_order = [B A C E D F]
  # chain_keys             = [0 1 2 3 4 5]
  consistent_chain_order = struct_asym_chain_ids.copy()
  consistent_chain_order[resolved_mask] = resolved_chain_ids
  return dict(zip(consistent_chain_order, range(len(struct_asym_chain_ids))))


def get_tables(
    cif: mmcif.Mmcif,
    fix_mse_residues: bool,
    fix_arginines: bool,
    fix_unknown_dna: bool,
    include_water: bool,
    include_other: bool,
    model_id: str,
) -> tuple[
    structure_tables.Chains, structure_tables.Residues, structure_tables.Atoms
]:
  """Returns chain, residue, and atom tables from a parsed mmcif.

  Args:
    cif: A parsed mmcif.Mmcif.
    fix_mse_residues: See from_mmcif.
    fix_arginines: See from_mmcif.
    fix_unknown_dna: See from_mmcif.
    include_water: See from_mmcif.
    include_other: See from_mmcif.
    model_id: A string defining which model ID to use. If set, only coordinates,
      b-factors and occupancies for the given model are returned. If empty,
      coordinates, b-factors and occupanciesall for models are returned with a
      leading dimension of num_models. Note that the model_id argument in
      from_mmcif is an integer and has slightly different use (see from_mmcif).
  """
  # Add any missing tables and columns we require for parsing.
  if cif_update := _generate_required_tables_if_missing(cif):
    cif = cif.copy_and_update(cif_update)

  # The wanted_chain_ids include chain IDs that are not necessarily present in
  # the _atom_site table, because we want to be able to read Structures with
  # chains consist of only unresolved residues.
  atom_site_all_models, wanted_chain_ids, layout = mmcif_utils.filter(
      cif,
      include_nucleotides=True,
      include_ligands=True,
      include_water=include_water,
      include_other=include_other,
      model_id=model_id,
  )
  atom_site_first_model = atom_site_all_models[0]

  # Get atom information from the _atom_site table.
  def _first_model_string_array(col: str) -> np.ndarray:
    return cif.get_array(col, dtype=object, gather=atom_site_first_model)

  def _requested_models_float_array(col: str) -> np.ndarray:
    if not model_id:
      # Return data for all models with a leading dimension of num_models.
      return cif.get_array(col, dtype=np.float32, gather=atom_site_all_models)
    else:
      # Return data only for the single requested model.
      return cif.get_array(col, dtype=np.float32, gather=atom_site_first_model)

  # These columns are the same for all models, fetch them just for the 1st one.
  label_comp_ids = _first_model_string_array('_atom_site.label_comp_id')
  label_asym_ids = _first_model_string_array('_atom_site.label_asym_id')
  label_seq_ids = _first_model_string_array('_atom_site.label_seq_id')
  label_atom_ids = _first_model_string_array('_atom_site.label_atom_id')
  if '_atom_site.auth_seq_id' in cif:
    auth_seq_ids = _first_model_string_array('_atom_site.auth_seq_id')
  else:
    auth_seq_ids = label_seq_ids  # auth_seq_id unset, fallback to label_seq_id.
  type_symbols = _first_model_string_array('_atom_site.type_symbol')
  pdbx_pdb_ins_codes = _first_model_string_array('_atom_site.pdbx_PDB_ins_code')

  # These columns are different for all models, fetch them as requested.
  atom_x = _requested_models_float_array('_atom_site.Cartn_x')
  atom_y = _requested_models_float_array('_atom_site.Cartn_y')
  atom_z = _requested_models_float_array('_atom_site.Cartn_z')
  atom_b_factor = _requested_models_float_array('_atom_site.B_iso_or_equiv')
  atom_occupancy = _requested_models_float_array('_atom_site.occupancy')

  # Make sure the scheme (residue) tables exist in case they are not present.
  if cif_update := _maybe_add_missing_scheme_tables(
      cif,
      res_starts=layout.residue_starts(),
      label_asym_ids=label_asym_ids,
      label_seq_ids=label_seq_ids,
      label_comp_ids=label_comp_ids,
      auth_seq_ids=auth_seq_ids,
      pdb_ins_codes=pdbx_pdb_ins_codes,
  ):
    cif = cif.copy_and_update(cif_update)

  # Fix common issues found in mmCIF files, like swapped arginine NH atoms.
  mmcif_utils.fix_residues(
      layout,
      comp_id=label_comp_ids,
      atom_id=label_atom_ids,
      atom_x=atom_x[0] if not model_id else atom_x,
      atom_y=atom_y[0] if not model_id else atom_y,
      atom_z=atom_z[0] if not model_id else atom_z,
      fix_arg=fix_arginines,
  )

  # Get keys for chains in the order they appear in _atom_site while also
  # dealing with empty chains.
  resolved_chain_ids = label_asym_ids[layout.chain_starts()]
  struct_asym_chain_ids = cif.get_array('_struct_asym.id', dtype=object)

  chain_key_by_chain_id = _get_chain_key_by_chain_id(
      resolved_chain_ids=resolved_chain_ids,
      struct_asym_chain_ids=struct_asym_chain_ids,
  )
  entity_id_by_chain_id = dict(
      zip(struct_asym_chain_ids, cif['_struct_asym.entity_id'], strict=True)
  )
  entity_description = cif.get(
      '_entity.pdbx_description', ['?'] * len(cif['_entity.id'])
  )
  entity_desc_by_entity_id = dict(
      zip(cif['_entity.id'], entity_description, strict=True)
  )
  chain_type_by_entity_id = mmcif.get_chain_type_by_entity_id(cif)
  auth_asym_id_by_chain_id = mmcif.get_internal_to_author_chain_id_map(cif)

  chain_res_builder = _ChainResBuilder(
      chain_key_by_chain_id=chain_key_by_chain_id,
      entity_id_by_chain_id=entity_id_by_chain_id,
      chain_type_by_entity_id=chain_type_by_entity_id,
      entity_desc_by_entity_id=entity_desc_by_entity_id,
      fix_mse_residues=fix_mse_residues,
      fix_unknown_dna=fix_unknown_dna,
  )

  # Collect data for polymer chain and residue tables. _pdbx_poly_seq_scheme is
  # guaranteed to be present thanks to _maybe_add_missing_scheme_tables.
  def _get_poly_seq_scheme_col(col: str) -> np.ndarray:
    return cif.get_array(key=f'_pdbx_poly_seq_scheme.{col}', dtype=object)

  poly_seq_asym_ids = _get_poly_seq_scheme_col('asym_id')
  poly_seq_pdb_seq_nums = _get_poly_seq_scheme_col('pdb_seq_num')
  poly_seq_seq_ids = _get_poly_seq_scheme_col('seq_id')
  poly_seq_mon_ids = _get_poly_seq_scheme_col('mon_id')
  poly_seq_pdb_strand_ids = _get_poly_seq_scheme_col('pdb_strand_id')
  poly_seq_pdb_ins_codes = _get_poly_seq_scheme_col('pdb_ins_code')
  string_array.remap(
      poly_seq_pdb_ins_codes, mapping=_INSERTION_CODE_REMAP, inplace=True
  )

  # We resolved alt-locs earlier for the atoms table. In cases of heterogeneous
  # residues (a residue with an alt-loc that is of different residue type), we
  # need to also do the same resolution in the residues table. Compute a mask
  # for the residues that were selected in the atoms table.
  poly_seq_mask = mmcif_utils.selected_polymer_residue_mask(
      layout=layout,
      atom_site_label_asym_ids=label_asym_ids[layout.residue_starts()],
      atom_site_label_seq_ids=label_seq_ids[layout.residue_starts()],
      atom_site_label_comp_ids=label_comp_ids[layout.residue_starts()],
      poly_seq_asym_ids=poly_seq_asym_ids,
      poly_seq_seq_ids=poly_seq_seq_ids,
      poly_seq_mon_ids=poly_seq_mon_ids,
  )

  if not include_other and poly_seq_mask:
    # Mask filtered-out residues so that they are not treated as missing.
    # Instead, we don't want them included in the chains/residues tables at all.
    keep_mask = string_array.isin(poly_seq_asym_ids, wanted_chain_ids)
    poly_seq_mask &= keep_mask

  chain_res_builder.add_residues(
      chain_ids=poly_seq_asym_ids[poly_seq_mask],
      chain_auth_asym_ids=poly_seq_pdb_strand_ids[poly_seq_mask],
      res_ids=poly_seq_seq_ids[poly_seq_mask].astype(np.int32),
      res_names=poly_seq_mon_ids[poly_seq_mask],
      res_auth_seq_ids=poly_seq_pdb_seq_nums[poly_seq_mask],
      res_ins_codes=poly_seq_pdb_ins_codes[poly_seq_mask],
  )

  # Collect data for ligand chain and residue tables. _pdbx_nonpoly_scheme
  # could be empty/unset if there are only branched ligands.
  def _get_nonpoly_scheme_col(col: str) -> np.ndarray:
    key = f'_pdbx_nonpoly_scheme.{col}'
    if f'_pdbx_nonpoly_scheme.{col}' in cif:
      return cif.get_array(key=key, dtype=object)
    else:
      return np.array([], dtype=object)

  nonpoly_asym_ids = _get_nonpoly_scheme_col('asym_id')
  nonpoly_auth_seq_ids = _get_nonpoly_scheme_col('pdb_seq_num')
  nonpoly_pdb_ins_codes = _get_nonpoly_scheme_col('pdb_ins_code')
  nonpoly_mon_ids = _get_nonpoly_scheme_col('mon_id')
  nonpoly_auth_asym_id = string_array.remap(
      nonpoly_asym_ids, mapping=auth_asym_id_by_chain_id, inplace=False
  )

  def _get_branch_scheme_col(col: str) -> np.ndarray:
    key = f'_pdbx_branch_scheme.{col}'
    if f'_pdbx_branch_scheme.{col}' in cif:
      return cif.get_array(key=key, dtype=object)
    else:
      return np.array([], dtype=object)

  branch_asym_ids = _get_branch_scheme_col('asym_id')
  branch_auth_seq_ids = _get_branch_scheme_col('pdb_seq_num')
  branch_pdb_ins_codes = _get_branch_scheme_col('pdb_ins_code')
  branch_mon_ids = _get_branch_scheme_col('mon_id')
  branch_auth_asym_id = string_array.remap(
      branch_asym_ids, mapping=auth_asym_id_by_chain_id, inplace=False
  )

  if branch_asym_ids.size > 0 and branch_pdb_ins_codes.size == 0:
    branch_pdb_ins_codes = np.array(['.'] * branch_asym_ids.size, dtype=object)

  # Compute the heterogeneous residue masks as above, this time for ligands.
  nonpoly_mask, branch_mask = mmcif_utils.selected_ligand_residue_mask(
      layout=layout,
      atom_site_label_asym_ids=label_asym_ids[layout.residue_starts()],
      atom_site_label_seq_ids=label_seq_ids[layout.residue_starts()],
      atom_site_auth_seq_ids=auth_seq_ids[layout.residue_starts()],
      atom_site_label_comp_ids=label_comp_ids[layout.residue_starts()],
      atom_site_pdbx_pdb_ins_codes=pdbx_pdb_ins_codes[layout.residue_starts()],
      nonpoly_asym_ids=nonpoly_asym_ids,
      nonpoly_auth_seq_ids=nonpoly_auth_seq_ids,
      nonpoly_pdb_ins_codes=nonpoly_pdb_ins_codes,
      nonpoly_mon_ids=nonpoly_mon_ids,
      branch_asym_ids=branch_asym_ids,
      branch_auth_seq_ids=branch_auth_seq_ids,
      branch_pdb_ins_codes=branch_pdb_ins_codes,
      branch_mon_ids=branch_mon_ids,
  )

  if not include_water:
    if nonpoly_mask:
      nonpoly_mask &= (nonpoly_mon_ids != 'HOH') & (nonpoly_mon_ids != 'DOD')
    if branch_mask:
      # Fix for bad mmCIFs that have water in the branch scheme table.
      branch_mask &= (branch_mon_ids != 'HOH') & (branch_mon_ids != 'DOD')

  string_array.remap(
      pdbx_pdb_ins_codes, mapping=_INSERTION_CODE_REMAP, inplace=True
  )
  string_array.remap(
      nonpoly_pdb_ins_codes, mapping=_INSERTION_CODE_REMAP, inplace=True
  )
  string_array.remap(
      branch_pdb_ins_codes, mapping=_INSERTION_CODE_REMAP, inplace=True
  )

  def _ligand_residue_ids(chain_ids: np.ndarray) -> np.ndarray:
    """Computes internal residue ID for ligand residues that don't have it."""

    # E.g. chain_ids=[A, A, A, B, C, C, D, D, D] -> [1, 2, 3, 1, 1, 2, 1, 2, 3].
    indices = np.arange(chain_ids.size, dtype=np.int32)
    return (indices + 1) - np.maximum.accumulate(
        indices * (chain_ids != np.roll(chain_ids, 1))
    )

  branch_residue_ids = _ligand_residue_ids(branch_asym_ids[branch_mask])
  nonpoly_residue_ids = _ligand_residue_ids(nonpoly_asym_ids[nonpoly_mask])

  chain_res_builder.add_residues(
      chain_ids=branch_asym_ids[branch_mask],
      chain_auth_asym_ids=branch_auth_asym_id[branch_mask],
      res_ids=branch_residue_ids,
      res_names=branch_mon_ids[branch_mask],
      res_auth_seq_ids=branch_auth_seq_ids[branch_mask],
      res_ins_codes=branch_pdb_ins_codes[branch_mask],
  )

  chain_res_builder.add_residues(
      chain_ids=nonpoly_asym_ids[nonpoly_mask],
      chain_auth_asym_ids=nonpoly_auth_asym_id[nonpoly_mask],
      res_ids=nonpoly_residue_ids,
      res_names=nonpoly_mon_ids[nonpoly_mask],
      res_auth_seq_ids=nonpoly_auth_seq_ids[nonpoly_mask],
      res_ins_codes=nonpoly_pdb_ins_codes[nonpoly_mask],
  )

  chains = chain_res_builder.make_chains_table()
  residues = chain_res_builder.make_residues_table()

  # Construct foreign residue keys for the atoms table.
  res_ends = np.array(layout.residues(), dtype=np.int32)
  res_starts = np.array(layout.residue_starts(), dtype=np.int32)
  res_lengths = res_ends - res_starts

  # Check just for HOH, DOD can be part e.g. of hydroxycysteine.
  if include_water:
    res_chain_types = chains.apply_array_to_column(
        column_name='type', arr=residues.chain_key
    )
    water_mask = res_chain_types != mmcif_names.WATER
    if 'HOH' in set(residues.name[water_mask]):
      raise ValueError('Bad mmCIF file: non-water entity has water molecules.')
  else:
    # Include resolved and unresolved residues.
    if 'HOH' in set(residues.name) | set(label_comp_ids[res_starts]):
      raise ValueError('Bad mmCIF file: non-water entity has water molecules.')

  atom_chain_key = string_array.remap(
      label_asym_ids, mapping=chain_res_builder.chain_key_by_chain_id
  ).astype(int)

  # If any of the residue lookups failed, the mmCIF is corrupted.
  try:
    atom_res_key_per_res = string_array.remap_multiple(
        (
            label_asym_ids[res_starts],
            auth_seq_ids[res_starts],
            label_comp_ids[res_starts],
            pdbx_pdb_ins_codes[res_starts],
        ),
        mapping=chain_res_builder.key_for_res,
    )
  except KeyError as e:
    raise ValueError(
        'Lookup for the following atom from the _atom_site table failed: '
        f'(label_asym_id, auth_seq_id, res_name, ins_code)={e}. This typically '
        'indicates that the _pdbx_poly_seq_scheme, _pdbx_nonpoly_scheme, or '
        '_pdbx_branch_scheme tables do not have data for all residues present '
        'in the _atom_site table. It could also be due to a known issue with '
        'a small number of multi-model mmCIFs.'
    ) from e

  # The residue ID will be shared for all atoms within that residue.
  atom_res_key = np.repeat(atom_res_key_per_res, repeats=res_lengths)

  if fix_mse_residues:
    met_residues_mask = (residues.name == 'MET')[atom_res_key]
    unfixed_mse_selenium_mask = met_residues_mask & (label_atom_ids == 'SE')
    label_atom_ids[unfixed_mse_selenium_mask] = 'SD'
    type_symbols[unfixed_mse_selenium_mask] = 'S'

  atoms = structure_tables.Atoms(
      key=atom_site_first_model,
      chain_key=atom_chain_key,
      res_key=atom_res_key,
      name=label_atom_ids,
      element=type_symbols,
      x=atom_x,
      y=atom_y,
      z=atom_z,
      b_factor=atom_b_factor,
      occupancy=atom_occupancy,
  )

  return chains, residues, atoms


def from_atom_arrays(
    *,
    res_id: np.ndarray,
    name: str = 'unset',
    release_date: datetime.date | None = None,
    resolution: float | None = None,
    structure_method: str | None = None,
    all_residues: Mapping[str, Sequence[tuple[str, int]]] | None = None,
    bioassembly_data: bioassemblies.BioassemblyData | None = None,
    chemical_components_data: (
        struc_chem_comps.ChemicalComponentsData | None
    ) = None,
    bond_table: structure_tables.Bonds | None = None,
    chain_id: np.ndarray | None = None,
    chain_type: np.ndarray | None = None,
    res_name: np.ndarray | None = None,
    atom_key: np.ndarray | None = None,
    atom_name: np.ndarray | None = None,
    atom_element: np.ndarray | None = None,
    atom_x: np.ndarray | None = None,
    atom_y: np.ndarray | None = None,
    atom_z: np.ndarray | None = None,
    atom_b_factor: np.ndarray | None = None,
    atom_occupancy: np.ndarray | None = None,
) -> structure.Structure:
  """Returns a Structure constructed from atom array level data.

  All fields except name and, res_id are optional, all array fields consist of a
  value for each atom in the structure - so residue and chain values should hold
  the same value for each atom in the chain or residue. Fields which are not
  defined are filled with default values.

  Validation is performed by the Structure constructor where possible - but
  author_naming scheme and all_residues must be checked in this function.

  It is not possible to construct structures with chains that do not contain
  any resolved residues using this function. If this is necessary, use the
  structure.Structure constructor directly.

  Args:
    res_id: Integer array of shape [num_atom]. The unique residue identifier for
      each residue. mmCIF field - _atom_site.label_seq_id.
    name: The name of the structure. E.g. a PDB ID.
    release_date: The release date of the structure as a `datetime.date`.
    resolution: The resolution of the structure in Angstroms.
    structure_method: The method used to solve this structure's coordinates.
    all_residues: An optional mapping from each chain ID (i.e. label_asym_id) to
      a sequence of (label_comp_id, label_seq_id) tuples, one per residue. This
      can contain residues that aren't present in the atom arrays. This is
      common in experimental data where some residues are not resolved but are
      known to be present.
    bioassembly_data: An optional instance of bioassembly.BioassemblyData. If
      present then a new Structure representing a specific bioassembly can be
      extracted using `Structure.generate_bioassembly(assembly_id)`.
    chemical_components_data: An optional instance of ChemicalComponentsData.
      Its content will be used for providing metadata about chemical components
      in this Structure instance. If not specified information will be retrieved
      from the standard chemical component dictionary (CCD, for more details see
      https://www.wwpdb.org/data/ccd).
    bond_table: A table representing manually-specified bonds. This corresponds
      to the _struct_conn table in an mmCIF. Atoms are identified by their key,
      as specified by the atom_key column. If this table is provided then the
      atom_key column must also be defined.
    chain_id: String array of shape [num_atom] of unique chain identifiers.
      mmCIF field - _atom_site.label_asym_id.
    chain_type: String array of shape [num_atom]. The molecular type of the
      current chain (e.g. polyribonucleotide). mmCIF field - _entity_poly.type
      OR _entity.type (for non-polymers).
    res_name: String array of shape [num_atom].. The name of each residue,
      typically a 3 letter string for polypeptides or 1-2 letter strings for
      polynucleotides. mmCIF field - _atom_site.label_comp_id.
    atom_key: A unique sorted integer array, used only by the bonds table to
      identify the atoms participating in each bond. If the bonds table is
      specified then this column must be non-None.
    atom_name: String array of shape [num_atom]. The name of each atom (e.g CA,
      O2', etc.). mmCIF field - _atom_site.label_atom_id.
    atom_element: String array of shape [num_atom]. The element type of each
      atom (e.g. C, O, N, etc.). mmCIF field - _atom_site.type_symbol.
    atom_x: Float array of shape [..., num_atom] of atom x coordinates. May have
      arbitrary leading dimensions, provided that these are consistent across
      all coordinate fields.
    atom_y: Float array of shape [..., num_atom] of atom y coordinates. May have
      arbitrary leading dimensions, provided that these are consistent across
      all coordinate fields.
    atom_z: Float array of shape [..., num_atom] of atom z coordinates. May have
      arbitrary leading dimensions, provided that these are consistent across
      all coordinate fields.
    atom_b_factor: Float array of shape [..., num_atom] or [num_atom] of atom
      b-factors or equivalent. If there are no extra leading dimensions then
      these values are assumed to apply to all coordinates for a given atom. If
      there are leading dimensions then these must match those used by the
      coordinate fields.
    atom_occupancy: Float array of shape [..., num_atom] or [num_atom] of atom
      occupancies or equivalent. If there are no extra leading dimensions then
      these values are assumed to apply to all coordinates for a given atom. If
      there are leading dimensions then these must match those used by the
      coordinate fields.
  """

  atoms, residues, chains = structure_tables.tables_from_atom_arrays(
      res_id=res_id,
      all_residues=all_residues,
      chain_id=chain_id,
      chain_type=chain_type,
      res_name=res_name,
      atom_key=atom_key,
      atom_name=atom_name,
      atom_element=atom_element,
      atom_x=atom_x,
      atom_y=atom_y,
      atom_z=atom_z,
      atom_b_factor=atom_b_factor,
      atom_occupancy=atom_occupancy,
  )

  return structure.Structure(
      name=name,
      release_date=release_date,
      resolution=resolution,
      structure_method=structure_method,
      bioassembly_data=bioassembly_data,
      chemical_components_data=chemical_components_data,
      atoms=atoms,
      chains=chains,
      residues=residues,
      bonds=bond_table or structure_tables.Bonds.make_empty(),
  )


def _guess_entity_type(
    chain_residues: Collection[str], atom_types: Collection[str]
) -> str:
  """Guess the entity type (polymer/non-polymer/water) based on residues/atoms.

  We treat both arguments as unordered collections since we care only whether
  all elements satisfy come conditions. The chain_residues can be either
  grouped by residue (length num_res), or it can be raw (length num_atoms).
  Atom type is unique for each atom in a residue, so don't group atom_types.

  Args:
    chain_residues: A sequence of full residue name (1-letter for DNA, 2-letters
      for RNA, 3 for protein). The _atom_site.label_comp_id column in mmCIF.
    atom_types: Atom type: ATOM or HETATM. The _atom_site.group_PDB column in
      mmCIF.

  Returns:
    One of polymer/non-polymer/water based on the following criteria:
    * If all atoms are HETATMs and all residues are water -> water.
    * If all atoms are HETATMs and not all residues are water -> non-polymer.
    * Otherwise -> polymer.
  """
  if not chain_residues or not atom_types:
    raise ValueError(
        f'chain_residues (len {len(chain_residues)}) and atom_types (len '
        f'{len(atom_types)}) must be both non-empty. Got: {chain_residues=} '
        f'and {atom_types=}'
    )

  if all(a == 'HETATM' for a in atom_types):
    if all(c in residue_names.WATER_TYPES for c in chain_residues):
      return mmcif_names.WATER
    return mmcif_names.NON_POLYMER_CHAIN
  return mmcif_names.POLYMER_CHAIN
