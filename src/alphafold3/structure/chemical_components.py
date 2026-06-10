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

"""Utilities for manipulating chemical components data."""

from collections.abc import Collection, Iterable, Mapping, Sequence
import dataclasses
import functools
from typing import Any, Self

from alphafold3.constants import chemical_components
from alphafold3.constants import residue_names
from alphafold3.structure import mmcif
import numpy as np
import rdkit.Chem as rd_chem


def _value_is_missing(value: Collection[Any] | str | None) -> bool:
  return value in ('.', '?', '', None)


def _to_optional_int(values: Sequence[str | None]) -> Sequence[int | None]:
  return [None if _value_is_missing(x) else int(x) for x in values]


def _to_optional_float(values: Sequence[str | None]) -> Sequence[float | None]:
  return [None if _value_is_missing(x) else float(x) for x in values]


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class ChemCompAtom:
  """Items of _chem_comp_atom category.

  See mmcif.wwpdb.org/dictionaries/mmcif_ma.dic/Categories/chem_comp_atom.html

  Attributes:
    type_symbol: _chem_comp_atom.type_symbol
    ordinal: _chem_comp_atom.pdbx_ordinal, can be optional
    charge: _chem_comp_atom.charge
    leaving_atom_flag: _chem_comp_atom.pdbx_leaving_atom_flag
    model_ideal_x: _chem_comp_atom.pdbx_model_Cartn_x_ideal
    model_ideal_y: _chem_comp_atom.pdbx_model_Cartn_y_ideal
    model_ideal_z: _chem_comp_atom.pdbx_model_Cartn_z_ideal
  """

  type_symbol: str
  ordinal: int | None = None
  charge: int | None = None
  leaving_atom_flag: str
  model_ideal_x: float | None = None
  model_ideal_y: float | None = None
  model_ideal_z: float | None = None


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class ChemCompBond:
  """Items of _chem_comp_bond category.

  See mmcif.wwpdb.org/dictionaries/mmcif_ma.dic/Categories/_chem_comp_bond.html

  Attributes:
    atom_id_1: _chem_comp_bond.atom_id_1
    atom_id_2: _chem_comp_bond.atom_id_2
    value_order: _chem_comp_bond.value_order
    aromatic_flag: _chem_comp_bond.aromatic_flag
    stereo_config: _chem_comp_bond.stereo_config
  """

  atom_id_1: str
  atom_id_2: str
  value_order: str
  aromatic_flag: str
  stereo_config: str


@dataclasses.dataclass(frozen=True)
class ChemCompEntry:
  """Items of _chem_comp category.

  For the full list of items and their semantics see
  http://mmcif.rcsb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/chem_comp.html
  """

  type: str
  name: str = '?'
  pdbx_synonyms: str = '?'
  formula: str = '?'
  formula_weight: str = '?'
  mon_nstd_flag: str = '?'
  pdbx_smiles: str | None = None
  # Mapping from _chem_comp_atom.atom_id to ChemCompAtom.
  chem_comp_atoms: Mapping[str, ChemCompAtom] | None = None
  chem_comp_bonds: list[ChemCompBond] | None = None

  def __post_init__(self):
    for field, value in vars(self).items():
      if not value and value is not None:
        raise ValueError(f"{field} value can't be an empty string.")

  def extends(self, other: Self) -> bool:
    """Checks whether this ChemCompEntry extends another one."""
    for field, value in vars(self).items():
      other_value = getattr(other, field)
      if _value_is_missing(other_value):
        continue
      if value != other_value:
        return False
    return True

  @property
  def rdkit_mol(self) -> rd_chem.Mol:
    """Returns an RDKit Mol, created via RDKit from entry SMILES string."""
    if not self.pdbx_smiles:
      raise ValueError('Cannot construct RDKit Mol with empty pdbx_smiles')
    return rd_chem.MolFromSmiles(self.pdbx_smiles)


_REQUIRED_MMCIF_COLUMNS = ('_chem_comp.id', '_chem_comp.type')


class MissingChemicalComponentsDataError(Exception):
  """Raised when chemical components data is missing from an mmCIF."""


@dataclasses.dataclass(frozen=True)
class ChemicalComponentsData:
  """Extra information for chemical components occurring in mmCIF.

  Fields:
    chem_comp: A mapping from _chem_comp.id to associated items in the
      chem_comp category.
  """

  chem_comp: Mapping[str, ChemCompEntry]

  @classmethod
  def from_mmcif(
      cls, cif: mmcif.Mmcif, fix_mse: bool, fix_unknown_dna: bool
  ) -> Self:
    """Constructs an instance of ChemicalComponentsData from an Mmcif object."""
    for col in _REQUIRED_MMCIF_COLUMNS:
      if col not in cif:
        raise MissingChemicalComponentsDataError(col)

    id_ = cif['_chem_comp.id']  # Guaranteed to be present.
    type_ = cif['_chem_comp.type']  # Guaranteed to be present.
    name = cif.get('_chem_comp.name', ['?'] * len(id_))
    synonyms = cif.get('_chem_comp.pdbx_synonyms', ['?'] * len(id_))
    formula = cif.get('_chem_comp.formula', ['?'] * len(id_))
    weight = cif.get('_chem_comp.formula_weight', ['?'] * len(id_))
    mon_nstd_flag = cif.get('_chem_comp.mon_nstd_flag', ['?'] * len(id_))
    smiles = cif.get('_chem_comp.pdbx_smiles', ['?'] * len(id_))
    smiles = [None if s == '?' else s for s in smiles]

    chem_comp_atoms_mapping = parse_atom_data(cif)
    chem_comp_bonds_mapping = parse_bond_data(cif)

    chem_comp = {}
    for component_name, *entry in zip(
        id_, type_, name, synonyms, formula, weight, mon_nstd_flag, smiles
    ):
      chem_comp_atoms = chem_comp_atoms_mapping.get(component_name)
      chem_comp_bonds = chem_comp_bonds_mapping.get(component_name)
      chem_comp[component_name] = ChemCompEntry(
          *entry,
          chem_comp_atoms=chem_comp_atoms,
          chem_comp_bonds=chem_comp_bonds,
      )

    if fix_mse and 'MSE' in chem_comp:
      if 'MET' not in chem_comp:
        chem_comp['MET'] = ChemCompEntry(
            type='L-PEPTIDE LINKING',
            name='METHIONINE',
            pdbx_synonyms='?',
            formula='C5 H11 N O2 S',
            formula_weight='149.211',
            mon_nstd_flag='y',
            pdbx_smiles=None,
        )

    if fix_unknown_dna and 'N' in chem_comp:
      # Do not delete 'N' as it may be needed for RNA in the system.
      if 'DN' not in chem_comp:
        chem_comp['DN'] = ChemCompEntry(
            type='DNA LINKING',
            name="UNKNOWN 2'-DEOXYNUCLEOTIDE",
            pdbx_synonyms='?',
            formula='C5 H11 O6 P',
            formula_weight='198.111',
            mon_nstd_flag='y',
            pdbx_smiles=None,
        )

    return ChemicalComponentsData(chem_comp)

  def to_mmcif_dict(self) -> Mapping[str, Sequence[str]]:
    """Returns chemical components data as a dict suitable for `mmcif.Mmcif`."""
    mmcif_dict = {}

    chem_comp_ids = []
    chem_comp_types = []
    chem_comp_names = []
    chem_comp_pdbx_synonyms = []
    chem_comp_formulas = []
    chem_comp_formula_weights = []
    chem_comp_mon_nstd_flags = []
    chem_comp_pdbx_smiles = []

    # _chem_comp_atom category
    chem_comp_atom_comp_ids = []
    chem_comp_atom_atom_ids = []
    chem_comp_atom_type_symbols = []
    chem_comp_atom_charges = []
    chem_comp_atom_leaving_atom_flags = []
    chem_comp_atom_model_ideal_xs = []
    chem_comp_atom_model_ideal_ys = []
    chem_comp_atom_model_ideal_zs = []
    chem_comp_atom_ordinals = []

    # _chem_comp_bond category
    chem_comp_bond_comp_ids = []
    chem_comp_bond_atom_id_1s = []
    chem_comp_bond_atom_id_2s = []
    chem_comp_bond_value_orders = []
    chem_comp_bond_aromatic_flags = []
    chem_comp_bond_stereo_configs = []

    for component_id in sorted(self.chem_comp):
      entry = self.chem_comp[component_id]
      chem_comp_ids.append(component_id)
      chem_comp_types.append(entry.type)
      chem_comp_names.append(entry.name)
      chem_comp_pdbx_synonyms.append(entry.pdbx_synonyms)
      chem_comp_formulas.append(entry.formula)
      chem_comp_formula_weights.append(entry.formula_weight)
      chem_comp_mon_nstd_flags.append(entry.mon_nstd_flag)
      chem_comp_pdbx_smiles.append(entry.pdbx_smiles)

      if entry.chem_comp_atoms:
        xs = []
        ys = []
        zs = []
        for atom_id, chem_comp_atom in entry.chem_comp_atoms.items():
          chem_comp_atom_comp_ids.append(component_id)
          chem_comp_atom_atom_ids.append(atom_id)
          chem_comp_atom_type_symbols.append(chem_comp_atom.type_symbol)
          chem_comp_atom_charges.append(
              str(chem_comp_atom.charge)
              if chem_comp_atom.charge is not None
              else '?'
          )
          chem_comp_atom_leaving_atom_flags.append(
              chem_comp_atom.leaving_atom_flag
          )
          xs.append(chem_comp_atom.model_ideal_x)
          ys.append(chem_comp_atom.model_ideal_y)
          zs.append(chem_comp_atom.model_ideal_z)
          chem_comp_atom_ordinals.append(
              str(chem_comp_atom.ordinal)
              if chem_comp_atom.ordinal is not None
              else '?'
          )

        def format_coords(coords: Sequence[float | None]) -> Sequence[str]:
          if any(c is None for c in coords):
            return [f'{c:.3f}' if c is not None else '?' for c in coords]
          else:
            return mmcif.format_float_array(
                np.array(coords, dtype=np.float64), num_decimal_places=3
            )

        chem_comp_atom_model_ideal_xs.extend(format_coords(xs))
        chem_comp_atom_model_ideal_ys.extend(format_coords(ys))
        chem_comp_atom_model_ideal_zs.extend(format_coords(zs))

      if entry.chem_comp_bonds:
        for bond in entry.chem_comp_bonds:
          chem_comp_bond_atom_id_1s.append(bond.atom_id_1)
          chem_comp_bond_atom_id_2s.append(bond.atom_id_2)
          chem_comp_bond_comp_ids.append(component_id)
          chem_comp_bond_value_orders.append(bond.value_order or '?')
          chem_comp_bond_aromatic_flags.append(bond.aromatic_flag or '?')
          chem_comp_bond_stereo_configs.append(bond.stereo_config or '?')

    if chem_comp_ids:
      mmcif_dict['_chem_comp.id'] = chem_comp_ids
      mmcif_dict['_chem_comp.type'] = chem_comp_types
      mmcif_dict['_chem_comp.name'] = chem_comp_names
      mmcif_dict['_chem_comp.pdbx_synonyms'] = chem_comp_pdbx_synonyms
      mmcif_dict['_chem_comp.formula'] = chem_comp_formulas
      mmcif_dict['_chem_comp.formula_weight'] = chem_comp_formula_weights
      mmcif_dict['_chem_comp.mon_nstd_flag'] = chem_comp_mon_nstd_flags

      if not all((v is None for v in chem_comp_pdbx_smiles)):
        mmcif_dict['_chem_comp.pdbx_smiles'] = [
            v or '?' for v in chem_comp_pdbx_smiles
        ]

    if chem_comp_atom_comp_ids:
      mmcif_dict['_chem_comp_atom.comp_id'] = chem_comp_atom_comp_ids
      mmcif_dict['_chem_comp_atom.atom_id'] = chem_comp_atom_atom_ids
      mmcif_dict['_chem_comp_atom.type_symbol'] = chem_comp_atom_type_symbols
      mmcif_dict['_chem_comp_atom.charge'] = chem_comp_atom_charges

      mmcif_dict['_chem_comp_atom.pdbx_leaving_atom_flag'] = (
          chem_comp_atom_leaving_atom_flags
      )

      mmcif_dict['_chem_comp_atom.pdbx_model_Cartn_x_ideal'] = (
          chem_comp_atom_model_ideal_xs
      )
      mmcif_dict['_chem_comp_atom.pdbx_model_Cartn_y_ideal'] = (
          chem_comp_atom_model_ideal_ys
      )
      mmcif_dict['_chem_comp_atom.pdbx_model_Cartn_z_ideal'] = (
          chem_comp_atom_model_ideal_zs
      )
      mmcif_dict['_chem_comp_atom.pdbx_ordinal'] = chem_comp_atom_ordinals

    if chem_comp_bond_comp_ids:
      mmcif_dict['_chem_comp_bond.comp_id'] = chem_comp_bond_comp_ids
      mmcif_dict['_chem_comp_bond.atom_id_1'] = chem_comp_bond_atom_id_1s
      mmcif_dict['_chem_comp_bond.atom_id_2'] = chem_comp_bond_atom_id_2s
      mmcif_dict['_chem_comp_bond.value_order'] = chem_comp_bond_value_orders
      mmcif_dict['_chem_comp_bond.pdbx_aromatic_flag'] = (
          chem_comp_bond_aromatic_flags
      )
      mmcif_dict['_chem_comp_bond.pdbx_stereo_config'] = (
          chem_comp_bond_stereo_configs
      )

    return mmcif_dict


def parse_atom_data(
    cif: mmcif.Mmcif,
) -> Mapping[str, Mapping[str, ChemCompAtom]]:
  """Parses _chem_comp_atom data from an Mmcif object."""

  atom_comps = cif.get('_chem_comp_atom.comp_id', [])
  num_atoms = len(atom_comps)
  if not num_atoms:
    return {}

  atom_ids = cif['_chem_comp_atom.atom_id']
  atom_types = cif['_chem_comp_atom.type_symbol']
  atom_leaving_atom_flags = cif.get(
      '_chem_comp_atom.pdbx_leaving_atom_flag', ['?'] * num_atoms
  )

  nones = [None] * num_atoms
  atom_charges = _to_optional_int(cif.get('_chem_comp_atom.charge', nones))
  atom_model_ideal_xs = _to_optional_float(
      cif.get('_chem_comp_atom.pdbx_model_Cartn_x_ideal', nones)
  )
  atom_model_ideal_ys = _to_optional_float(
      cif.get('_chem_comp_atom.pdbx_model_Cartn_y_ideal', nones)
  )
  atom_model_ideal_zs = _to_optional_float(
      cif.get('_chem_comp_atom.pdbx_model_Cartn_z_ideal', nones)
  )
  ordinal = _to_optional_int(cif.get('_chem_comp_atom.pdbx_ordinal', nones))

  chem_comp_atoms = {}
  for i, atom_comp in enumerate(atom_comps):
    comp_atom_vals = chem_comp_atoms.setdefault(atom_comp, {})
    comp_atom_vals[atom_ids[i]] = ChemCompAtom(
        type_symbol=atom_types[i],
        ordinal=ordinal[i],
        charge=atom_charges[i],
        leaving_atom_flag=atom_leaving_atom_flags[i],
        model_ideal_x=atom_model_ideal_xs[i],
        model_ideal_y=atom_model_ideal_ys[i],
        model_ideal_z=atom_model_ideal_zs[i],
    )
  return chem_comp_atoms


def parse_bond_data(
    cif: mmcif.Mmcif,
) -> Mapping[str, list[ChemCompBond]]:
  """Parses _chem_comp_bond data from an Mmcif object."""
  bond_comps = cif.get('_chem_comp_bond.comp_id', [])
  len_bonds = len(bond_comps)
  if not len_bonds:
    return {}

  bond_atom_ids_1 = cif['_chem_comp_bond.atom_id_1']
  bond_atom_ids_2 = cif['_chem_comp_bond.atom_id_2']

  unknowns = ['?'] * len_bonds
  bond_value_orders = cif.get('_chem_comp_bond.value_order', unknowns)
  bond_aromatic_flags = cif.get('_chem_comp_bond.pdbx_aromatic_flag', unknowns)
  bond_stereo_configs = cif.get('_chem_comp_bond.pdbx_stereo_config', unknowns)

  chem_comp_bonds = {}
  for i, bond_comp in enumerate(bond_comps):
    chem_comp_bonds.setdefault(bond_comp, []).append(
        ChemCompBond(
            atom_id_1=bond_atom_ids_1[i],
            atom_id_2=bond_atom_ids_2[i],
            value_order=bond_value_orders[i],
            aromatic_flag=bond_aromatic_flags[i],
            stereo_config=bond_stereo_configs[i],
        )
    )
  return chem_comp_bonds


def get_data_for_ccd_components(
    ccd: chemical_components.Ccd,
    chemical_component_ids: Iterable[str],
    populate_pdbx_smiles: bool = False,
) -> ChemicalComponentsData:
  """Returns `ChemicalComponentsData` for chemical components known by PDB."""
  chem_comp = {}
  for chemical_component_id in chemical_component_ids:
    chem_data = chemical_components.component_name_to_info(
        ccd=ccd, res_name=chemical_component_id
    )
    if not chem_data:
      continue

    chem_comp[chemical_component_id] = ChemCompEntry(
        type=chem_data.type,
        name=chem_data.name,
        pdbx_synonyms=chem_data.pdbx_synonyms,
        formula=chem_data.formula,
        formula_weight=chem_data.formula_weight,
        mon_nstd_flag=chem_data.mon_nstd_flag,
        pdbx_smiles=(
            chem_data.pdbx_smiles or None if populate_pdbx_smiles else None
        ),
    )
  return ChemicalComponentsData(chem_comp=chem_comp)


def populate_missing_ccd_data(
    ccd: chemical_components.Ccd,
    chemical_components_data: ChemicalComponentsData,
    chemical_component_ids: Iterable[str] | None = None,
    populate_pdbx_smiles: bool = False,
) -> ChemicalComponentsData:
  """Populates missing data for the chemical components from CCD.

  Args:
    ccd: The chemical components database.
    chemical_components_data: ChemicalComponentsData to populate missing values
      for. This function doesn't modify the object, extended version is provided
      as a return value.
    chemical_component_ids: chemical components to populate missing values for.
      If not specified, the function will consider all chemical components which
      are already present in `chemical_components_data`.
    populate_pdbx_smiles: whether to populate `pdbx_smiles` field using SMILES
      descriptors from _pdbx_chem_comp_descriptor CCD table. If CCD provides
      multiple SMILES strings, any of them could be used.

  Returns:
    New instance of ChemicalComponentsData without missing values for CCD
    entries.
  """
  if chemical_component_ids is None:
    chemical_component_ids = chemical_components_data.chem_comp.keys()

  ccd_data = get_data_for_ccd_components(
      ccd, chemical_component_ids, populate_pdbx_smiles
  )
  chem_comp = dict(chemical_components_data.chem_comp)
  for component_id, ccd_entry in ccd_data.chem_comp.items():
    if component_id not in chem_comp:
      chem_comp[component_id] = ccd_entry
    else:
      already_specified_fields = {
          field: value
          for field, value in vars(chem_comp[component_id]).items()
          if not _value_is_missing(value)
      }
      chem_comp[component_id] = ChemCompEntry(
          **{**vars(ccd_entry), **already_specified_fields}
      )
  return ChemicalComponentsData(chem_comp=chem_comp)


def get_all_atoms_in_entry(
    ccd: chemical_components.Ccd, res_name: str
) -> Mapping[str, Sequence[str]]:
  """Get all possible atoms and bonds for this residue in a standard order.

  Args:
    ccd: The chemical components dictionary.
    res_name: Full CCD name.

  Returns:
    A dictionary table of the atoms and bonds for this residue in this residue
    type.
  """
  # The CCD version of 'UNK' is weird. It has a CB and a CG atom. We just want
  # the minimal amino-acid here which is GLY.
  if res_name == 'UNK':
    res_name = 'GLY'
  ccd_data = ccd.get(res_name)
  if not ccd_data:
    raise ValueError(f'Unknown residue type {res_name}')

  keys = (
      '_chem_comp_atom.atom_id',
      '_chem_comp_atom.type_symbol',
      '_chem_comp_bond.atom_id_1',
      '_chem_comp_bond.atom_id_2',
  )

  # Add terminal hydrogens for protonation of the N-terminal
  if res_name == 'PRO':
    res_atoms = {key: [*ccd_data.get(key, [])] for key in keys}
    res_atoms['_chem_comp_atom.atom_id'].extend(['H2', 'H3'])
    res_atoms['_chem_comp_atom.type_symbol'].extend(['H', 'H'])
    res_atoms['_chem_comp_bond.atom_id_1'].extend(['N', 'N'])
    res_atoms['_chem_comp_bond.atom_id_2'].extend(['H2', 'H3'])
  elif res_name in residue_names.PROTEIN_TYPES_WITH_UNKNOWN:
    res_atoms = {key: [*ccd_data.get(key, [])] for key in keys}
    res_atoms['_chem_comp_atom.atom_id'].append('H3')
    res_atoms['_chem_comp_atom.type_symbol'].append('H')
    res_atoms['_chem_comp_bond.atom_id_1'].append('N')
    res_atoms['_chem_comp_bond.atom_id_2'].append('H3')
  else:
    res_atoms = {key: ccd_data.get(key, []) for key in keys}

  return res_atoms


@functools.lru_cache(maxsize=128)
def get_res_atom_names(ccd: chemical_components.Ccd, res_name: str) -> set[str]:
  """Gets the names of the atoms in a given CCD residue."""
  atoms = get_all_atoms_in_entry(ccd, res_name)['_chem_comp_atom.atom_id']
  return set(atoms)
