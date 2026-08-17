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

"""List of atom types with reverse look-up."""

from collections.abc import Mapping, Sequence, Set
import itertools
import sys
from typing import Final
from alphafold3.constants import residue_names

# Note:
# `sys.intern` places the values in the Python internal db for fast lookup.

# 37 common residue atoms.
N = sys.intern('N')
CA = sys.intern('CA')
C = sys.intern('C')
CB = sys.intern('CB')
O = sys.intern('O')
CG = sys.intern('CG')
CG1 = sys.intern('CG1')
CG2 = sys.intern('CG2')
OG = sys.intern('OG')
OG1 = sys.intern('OG1')
SG = sys.intern('SG')
CD = sys.intern('CD')
CD1 = sys.intern('CD1')
CD2 = sys.intern('CD2')
ND1 = sys.intern('ND1')
ND2 = sys.intern('ND2')
OD1 = sys.intern('OD1')
OD2 = sys.intern('OD2')
SD = sys.intern('SD')
CE = sys.intern('CE')
CE1 = sys.intern('CE1')
CE2 = sys.intern('CE2')
CE3 = sys.intern('CE3')
NE = sys.intern('NE')
NE1 = sys.intern('NE1')
NE2 = sys.intern('NE2')
OE1 = sys.intern('OE1')
OE2 = sys.intern('OE2')
CH2 = sys.intern('CH2')
NH1 = sys.intern('NH1')
NH2 = sys.intern('NH2')
OH = sys.intern('OH')
CZ = sys.intern('CZ')
CZ2 = sys.intern('CZ2')
CZ3 = sys.intern('CZ3')
NZ = sys.intern('NZ')
OXT = sys.intern('OXT')

# 29 common nucleic acid atoms.
C1PRIME = sys.intern("C1'")
C2 = sys.intern('C2')
C2PRIME = sys.intern("C2'")
C3PRIME = sys.intern("C3'")
C4 = sys.intern('C4')
C4PRIME = sys.intern("C4'")
C5 = sys.intern('C5')
C5PRIME = sys.intern("C5'")
C6 = sys.intern('C6')
C7 = sys.intern('C7')
C8 = sys.intern('C8')
N1 = sys.intern('N1')
N2 = sys.intern('N2')
N3 = sys.intern('N3')
N4 = sys.intern('N4')
N6 = sys.intern('N6')
N7 = sys.intern('N7')
N9 = sys.intern('N9')
O2 = sys.intern('O2')
O2PRIME = sys.intern("O2'")
O3PRIME = sys.intern("O3'")
O4 = sys.intern('O4')
O4PRIME = sys.intern("O4'")
O5PRIME = sys.intern("O5'")
O6 = sys.intern('O6')
OP1 = sys.intern('OP1')
OP2 = sys.intern('OP2')
OP3 = sys.intern('OP3')
P = sys.intern('P')

# A list of atoms (excluding hydrogen) for each AA type. PDB naming convention.
RESIDUE_ATOMS: Mapping[str, tuple[str, ...]] = {
    residue_names.ALA: (C, CA, CB, N, O),
    residue_names.ARG: (C, CA, CB, CG, CD, CZ, N, NE, O, NH1, NH2),
    residue_names.ASN: (C, CA, CB, CG, N, ND2, O, OD1),
    residue_names.ASP: (C, CA, CB, CG, N, O, OD1, OD2),
    residue_names.CYS: (C, CA, CB, N, O, SG),
    residue_names.GLN: (C, CA, CB, CG, CD, N, NE2, O, OE1),
    residue_names.GLU: (C, CA, CB, CG, CD, N, O, OE1, OE2),
    residue_names.GLY: (C, CA, N, O),
    residue_names.HIS: (C, CA, CB, CG, CD2, CE1, N, ND1, NE2, O),
    residue_names.ILE: (C, CA, CB, CG1, CG2, CD1, N, O),
    residue_names.LEU: (C, CA, CB, CG, CD1, CD2, N, O),
    residue_names.LYS: (C, CA, CB, CG, CD, CE, N, NZ, O),
    residue_names.MET: (C, CA, CB, CG, CE, N, O, SD),
    residue_names.PHE: (C, CA, CB, CG, CD1, CD2, CE1, CE2, CZ, N, O),
    residue_names.PRO: (C, CA, CB, CG, CD, N, O),
    residue_names.SER: (C, CA, CB, N, O, OG),
    residue_names.THR: (C, CA, CB, CG2, N, O, OG1),
    residue_names.TRP:
        (C, CA, CB, CG, CD1, CD2, CE2, CE3, CZ2, CZ3, CH2, N, NE1, O),
    residue_names.TYR: (C, CA, CB, CG, CD1, CD2, CE1, CE2, CZ, N, O, OH),
    residue_names.VAL: (C, CA, CB, CG1, CG2, N, O),
}  # pyformat: disable

# Used to identify backbone for alignment and distance calculation for sterics.
PROTEIN_BACKBONE_ATOMS: tuple[str, ...] = (N, CA, C)

# Naming swaps for ambiguous atom names. Due to symmetries in the amino acids
# the naming of atoms is ambiguous in 4 of the 20 amino acids. (The LDDT paper
# lists 7 amino acids as ambiguous, but the naming ambiguities in LEU, VAL and
# ARG can be resolved by using the 3D constellations of the 'ambiguous' atoms
# and their neighbours)
AMBIGUOUS_ATOM_NAMES: Mapping[str, Mapping[str, str]] = {
    residue_names.ASP: {OD1: OD2},
    residue_names.GLU: {OE1: OE2},
    residue_names.PHE: {CD1: CD2, CE1: CE2},
    residue_names.TYR: {CD1: CD2, CE1: CE2},
}

# Used when we need to store atom data in a format that requires fixed atom data
# size for every protein residue (e.g. a numpy array).
ATOM37: tuple[str, ...] = (
    N, CA, C, CB, O, CG, CG1, CG2, OG, OG1, SG, CD, CD1, CD2, ND1, ND2, OD1,
    OD2, SD, CE, CE1, CE2, CE3, NE, NE1, NE2, OE1, OE2, CH2, NH1, NH2, OH, CZ,
    CZ2, CZ3, NZ, OXT)  # pyformat: disable
ATOM37_ORDER: Mapping[str, int] = {name: i for i, name in enumerate(ATOM37)}
ATOM37_NUM: Final[int] = len(ATOM37)  # := 37.

# Used when we need to store protein atom data in a format that requires fixed
# atom data size for any residue but takes less space than ATOM37 by having 14
# fields, which is sufficient for storing atoms of all protein residues (e.g. a
# numpy array).
ATOM14: Mapping[str, tuple[str, ...]] = {
    residue_names.ALA: (N, CA, C, O, CB),
    residue_names.ARG: (N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2),
    residue_names.ASN: (N, CA, C, O, CB, CG, OD1, ND2),
    residue_names.ASP: (N, CA, C, O, CB, CG, OD1, OD2),
    residue_names.CYS: (N, CA, C, O, CB, SG),
    residue_names.GLN: (N, CA, C, O, CB, CG, CD, OE1, NE2),
    residue_names.GLU: (N, CA, C, O, CB, CG, CD, OE1, OE2),
    residue_names.GLY: (N, CA, C, O),
    residue_names.HIS: (N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2),
    residue_names.ILE: (N, CA, C, O, CB, CG1, CG2, CD1),
    residue_names.LEU: (N, CA, C, O, CB, CG, CD1, CD2),
    residue_names.LYS: (N, CA, C, O, CB, CG, CD, CE, NZ),
    residue_names.MET: (N, CA, C, O, CB, CG, SD, CE),
    residue_names.PHE: (N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ),
    residue_names.PRO: (N, CA, C, O, CB, CG, CD),
    residue_names.SER: (N, CA, C, O, CB, OG),
    residue_names.THR: (N, CA, C, O, CB, OG1, CG2),
    residue_names.TRP:
        (N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2),
    residue_names.TYR: (N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH),
    residue_names.VAL: (N, CA, C, O, CB, CG1, CG2),
    residue_names.UNK: (),
}  # pyformat: disable

# A compact atom encoding with 14 columns, padded with '' in empty slots.
ATOM14_PADDED: Mapping[str, Sequence[str]] = {
    k: [v for _, v in itertools.zip_longest(range(14), values, fillvalue='')]
    for k, values in ATOM14.items()
}

ATOM14_ORDER: Mapping[str, Mapping[str, int]] = {
    k: {name: i for i, name in enumerate(v)} for k, v in ATOM14.items()
}
ATOM14_NUM: Final[int] = max(len(v) for v in ATOM14.values())

# Used when we need to store protein and nucleic atom library.
DENSE_ATOM: Mapping[str, tuple[str, ...]] = {
    # Protein.
    residue_names.ALA: (N, CA, C, O, CB),
    residue_names.ARG: (N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2),
    residue_names.ASN: (N, CA, C, O, CB, CG, OD1, ND2),
    residue_names.ASP: (N, CA, C, O, CB, CG, OD1, OD2),
    residue_names.CYS: (N, CA, C, O, CB, SG),
    residue_names.GLN: (N, CA, C, O, CB, CG, CD, OE1, NE2),
    residue_names.GLU: (N, CA, C, O, CB, CG, CD, OE1, OE2),
    residue_names.GLY: (N, CA, C, O),
    residue_names.HIS: (N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2),
    residue_names.ILE: (N, CA, C, O, CB, CG1, CG2, CD1),
    residue_names.LEU: (N, CA, C, O, CB, CG, CD1, CD2),
    residue_names.LYS: (N, CA, C, O, CB, CG, CD, CE, NZ),
    residue_names.MET: (N, CA, C, O, CB, CG, SD, CE),
    residue_names.PHE: (N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ),
    residue_names.PRO: (N, CA, C, O, CB, CG, CD),
    residue_names.SER: (N, CA, C, O, CB, OG),
    residue_names.THR: (N, CA, C, O, CB, OG1, CG2),
    residue_names.TRP:
        (N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2),
    residue_names.TYR: (N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH),
    residue_names.VAL: (N, CA, C, O, CB, CG1, CG2),
    residue_names.UNK: (),
    # RNA.
    residue_names.A:
        (OP3, P, OP1, OP2, O5PRIME, C5PRIME, C4PRIME, O4PRIME, C3PRIME, O3PRIME,
         C2PRIME, O2PRIME, C1PRIME, N9, C8, N7, C5, C6, N6, N1, C2, N3, C4),
    residue_names.C:
        (OP3, P, OP1, OP2, O5PRIME, C5PRIME, C4PRIME, O4PRIME, C3PRIME, O3PRIME,
         C2PRIME, O2PRIME, C1PRIME, N1, C2, O2, N3, C4, N4, C5, C6),
    residue_names.G:
        (OP3, P, OP1, OP2, O5PRIME, C5PRIME, C4PRIME, O4PRIME, C3PRIME, O3PRIME,
         C2PRIME, O2PRIME, C1PRIME, N9, C8, N7, C5, C6, O6, N1, C2, N2, N3, C4),
    residue_names.U:
        (OP3, P, OP1, OP2, O5PRIME, C5PRIME, C4PRIME, O4PRIME, C3PRIME, O3PRIME,
         C2PRIME, O2PRIME, C1PRIME, N1, C2, O2, N3, C4, O4, C5, C6),
    residue_names.UNK_RNA: (),
    # DNA.
    residue_names.DA:
        (OP3, P, OP1, OP2, O5PRIME, C5PRIME, C4PRIME, O4PRIME, C3PRIME, O3PRIME,
         C2PRIME, C1PRIME, N9, C8, N7, C5, C6, N6, N1, C2, N3, C4),
    residue_names.DC:
        (OP3, P, OP1, OP2, O5PRIME, C5PRIME, C4PRIME, O4PRIME, C3PRIME, O3PRIME,
         C2PRIME, C1PRIME, N1, C2, O2, N3, C4, N4, C5, C6),
    residue_names.DG:
        (OP3, P, OP1, OP2, O5PRIME, C5PRIME, C4PRIME, O4PRIME, C3PRIME, O3PRIME,
         C2PRIME, C1PRIME, N9, C8, N7, C5, C6, O6, N1, C2, N2, N3, C4),
    residue_names.DT:
        (OP3, P, OP1, OP2, O5PRIME, C5PRIME, C4PRIME, O4PRIME, C3PRIME, O3PRIME,
         C2PRIME, C1PRIME, N1, C2, O2, N3, C4, O4, C5, C7, C6),
    # Unknown nucleic.
    residue_names.UNK_DNA: (),
}  # pyformat: disable

DENSE_ATOM_ORDER: Mapping[str, Mapping[str, int]] = {
    k: {name: i for i, name in enumerate(v)} for k, v in DENSE_ATOM.items()
}
DENSE_ATOM_NUM: Final[int] = max(len(v) for v in DENSE_ATOM.values())

# Used when we need to store atom data in a format that requires fixed atom data
# size for every nucleic molecule (e.g. a numpy array).
ATOM29: tuple[str, ...] = (
    "C1'", 'C2', "C2'", "C3'", 'C4', "C4'", 'C5', "C5'", 'C6', 'C7', 'C8', 'N1',
    'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'OP3', 'O2', "O2'", "O3'", 'O4', "O4'",
    "O5'", 'O6', 'OP1', 'OP2', 'P')  # pyformat: disable
ATOM29_ORDER: Mapping[str, int] = {
    atom_type: i for i, atom_type in enumerate(ATOM29)
}
ATOM29_NUM: Final[int] = len(ATOM29)  # := 29

# Hydrogens that exist depending on the protonation state of the residue.
# Extracted from third_party/py/openmm/app/data/hydrogens.xml
PROTONATION_HYDROGENS: Mapping[str, Set[str]] = {
    'ASP': {'HD2'},
    'CYS': {'HG'},
    'GLU': {'HE2'},
    'HIS': {'HD1', 'HE2'},
    'LYS': {'HZ3'},
}
