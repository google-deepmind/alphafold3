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

from collections.abc import Sequence

import numpy as np

from alphafold3.cpp import cif_dict
from alphafold3.structure.python import mmcif_layout


def filter(
    mmcif: cif_dict.CifDict,
    include_nucleotides: bool,
    include_ligands: bool = ...,
    include_water: bool = ...,
    include_other: bool = ...,
    model_id: str = ...,
) -> tuple[np.ndarray[int], set[str], mmcif_layout.MmcifLayout]: ...


def fix_residues(
    layout: mmcif_layout.MmcifLayout,
    comp_id: Sequence[str],
    atom_id: Sequence[str],
    atom_x: Sequence[float],
    atom_y: Sequence[float],
    atom_z: Sequence[float],
    fix_arg: bool = ...,
) -> None: ...


def map_residue_keys(
    src_arrays: Sequence[np.ndarray],
    query_arrays: Sequence[np.ndarray],
) -> np.ndarray: ...
def read_layout(
    mmcif: cif_dict.CifDict, model_id: str = ...
) -> mmcif_layout.MmcifLayout: ...


def selected_ligand_residue_mask(
    layout: mmcif_layout.MmcifLayout,
    atom_site_label_asym_ids: Sequence[str] | np.ndarray,
    atom_site_label_seq_ids: Sequence[str] | np.ndarray,
    atom_site_auth_seq_ids: Sequence[str] | np.ndarray,
    atom_site_label_comp_ids: Sequence[str] | np.ndarray,
    atom_site_pdbx_pdb_ins_codes: Sequence[str] | np.ndarray,
    nonpoly_asym_ids: Sequence[str] | np.ndarray,
    nonpoly_auth_seq_ids: Sequence[str] | np.ndarray,
    nonpoly_pdb_ins_codes: Sequence[str] | np.ndarray,
    nonpoly_mon_ids: Sequence[str] | np.ndarray,
    branch_asym_ids: Sequence[str] | np.ndarray,
    branch_auth_seq_ids: Sequence[str] | np.ndarray,
    branch_pdb_ins_codes: Sequence[str] | np.ndarray,
    branch_mon_ids: Sequence[str] | np.ndarray,
) -> tuple[np.ndarray, np.ndarray]: ...


def selected_polymer_residue_mask(
    layout: mmcif_layout.MmcifLayout,
    atom_site_label_asym_ids: Sequence[str] | np.ndarray,
    atom_site_label_seq_ids: Sequence[str] | np.ndarray,
    atom_site_label_comp_ids: Sequence[str] | np.ndarray,
    poly_seq_asym_ids: Sequence[str] | np.ndarray,
    poly_seq_seq_ids: Sequence[str] | np.ndarray,
    poly_seq_mon_ids: Sequence[str] | np.ndarray,
) -> np.ndarray: ...
