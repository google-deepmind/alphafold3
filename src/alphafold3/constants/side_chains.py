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

"""Constants associated with side chains."""

from collections.abc import Mapping, Sequence
import itertools

# Format: The list for each AA type contains chi1, chi2, chi3, chi4 in
# this order (or a relevant subset from chi1 onwards). ALA and GLY don't have
# chi angles so their chi angle lists are empty.
CHI_ANGLES_ATOMS: Mapping[str, Sequence[tuple[str, ...]]] = {
    'ALA': [],
    # Chi5 in arginine is always 0 +- 5 degrees, so ignore it.
    'ARG': [
        ('N', 'CA', 'CB', 'CG'),
        ('CA', 'CB', 'CG', 'CD'),
        ('CB', 'CG', 'CD', 'NE'),
        ('CG', 'CD', 'NE', 'CZ'),
    ],
    'ASN': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'OD1')],
    'ASP': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'OD1')],
    'CYS': [('N', 'CA', 'CB', 'SG')],
    'GLN': [
        ('N', 'CA', 'CB', 'CG'),
        ('CA', 'CB', 'CG', 'CD'),
        ('CB', 'CG', 'CD', 'OE1'),
    ],
    'GLU': [
        ('N', 'CA', 'CB', 'CG'),
        ('CA', 'CB', 'CG', 'CD'),
        ('CB', 'CG', 'CD', 'OE1'),
    ],
    'GLY': [],
    'HIS': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'ND1')],
    'ILE': [('N', 'CA', 'CB', 'CG1'), ('CA', 'CB', 'CG1', 'CD1')],
    'LEU': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'LYS': [
        ('N', 'CA', 'CB', 'CG'),
        ('CA', 'CB', 'CG', 'CD'),
        ('CB', 'CG', 'CD', 'CE'),
        ('CG', 'CD', 'CE', 'NZ'),
    ],
    'MET': [
        ('N', 'CA', 'CB', 'CG'),
        ('CA', 'CB', 'CG', 'SD'),
        ('CB', 'CG', 'SD', 'CE'),
    ],
    'PHE': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'PRO': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')],
    'SER': [('N', 'CA', 'CB', 'OG')],
    'THR': [('N', 'CA', 'CB', 'OG1')],
    'TRP': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'TYR': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'VAL': [('N', 'CA', 'CB', 'CG1')],
}

_CHI_GROUPS_FOR_ATOM: dict[tuple[str, str], list[tuple[int, int]]] = {}
for res_name, chi_angle_atoms_for_res in CHI_ANGLES_ATOMS.items():
  for chi_group_i, chi_group in enumerate(chi_angle_atoms_for_res):
    for atom_i, atom in enumerate(chi_group):
      _CHI_GROUPS_FOR_ATOM.setdefault((res_name, atom), []).append(
          (chi_group_i, atom_i)
      )

# Mapping from (residue_name, atom_name) pairs to the atom's chi group index
# and atom index within that group.
CHI_GROUPS_FOR_ATOM: Mapping[tuple[str, str], Sequence[tuple[int, int]]] = (
    _CHI_GROUPS_FOR_ATOM
)

MAX_NUM_CHI_ANGLES: int = 4
ATOMS_PER_CHI_ANGLE: int = 4

# A list of atoms for each AA type that are involved in chi angle calculations.
CHI_ATOM_SETS: Mapping[str, set[str]] = {
    residue_name: set(itertools.chain(*atoms))
    for residue_name, atoms in CHI_ANGLES_ATOMS.items()
}

# If chi angles given in fixed-length array, this matrix determines how to mask
# them for each AA type. The order is as per restype_order (see below).
CHI_ANGLES_MASK: Sequence[Sequence[float]] = (
    (0.0, 0.0, 0.0, 0.0),  # ALA
    (1.0, 1.0, 1.0, 1.0),  # ARG
    (1.0, 1.0, 0.0, 0.0),  # ASN
    (1.0, 1.0, 0.0, 0.0),  # ASP
    (1.0, 0.0, 0.0, 0.0),  # CYS
    (1.0, 1.0, 1.0, 0.0),  # GLN
    (1.0, 1.0, 1.0, 0.0),  # GLU
    (0.0, 0.0, 0.0, 0.0),  # GLY
    (1.0, 1.0, 0.0, 0.0),  # HIS
    (1.0, 1.0, 0.0, 0.0),  # ILE
    (1.0, 1.0, 0.0, 0.0),  # LEU
    (1.0, 1.0, 1.0, 1.0),  # LYS
    (1.0, 1.0, 1.0, 0.0),  # MET
    (1.0, 1.0, 0.0, 0.0),  # PHE
    (1.0, 1.0, 0.0, 0.0),  # PRO
    (1.0, 0.0, 0.0, 0.0),  # SER
    (1.0, 0.0, 0.0, 0.0),  # THR
    (1.0, 1.0, 0.0, 0.0),  # TRP
    (1.0, 1.0, 0.0, 0.0),  # TYR
    (1.0, 0.0, 0.0, 0.0),  # VAL
)
