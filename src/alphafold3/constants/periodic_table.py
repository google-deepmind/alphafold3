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

"""Periodic table of elements."""

from collections.abc import Mapping, Sequence
import dataclasses
from typing import Final

import numpy as np


@dataclasses.dataclass(frozen=True, kw_only=True)
class Element:
  name: str
  number: int
  symbol: str
  weight: float


# Weights taken from rdkit/Code/GraphMol/atomic_data.cpp for compatibility.
# pylint: disable=invalid-name

# X is an unknown element that can be present in the CCD,
# https://www.rcsb.org/ligand/UNX.
X: Final[Element] = Element(name='Unknown', number=0, symbol='X', weight=0.0)
H: Final[Element] = Element(name='Hydrogen', number=1, symbol='H', weight=1.008)
He: Final[Element] = Element(name='Helium', number=2, symbol='He', weight=4.003)
Li: Final[Element] = Element(
    name='Lithium', number=3, symbol='Li', weight=6.941
)
Be: Final[Element] = Element(
    name='Beryllium', number=4, symbol='Be', weight=9.012
)
B: Final[Element] = Element(name='Boron', number=5, symbol='B', weight=10.812)
C: Final[Element] = Element(name='Carbon', number=6, symbol='C', weight=12.011)
N: Final[Element] = Element(
    name='Nitrogen', number=7, symbol='N', weight=14.007
)
O: Final[Element] = Element(name='Oxygen', number=8, symbol='O', weight=15.999)
F: Final[Element] = Element(
    name='Fluorine', number=9, symbol='F', weight=18.998
)
Ne: Final[Element] = Element(name='Neon', number=10, symbol='Ne', weight=20.18)
Na: Final[Element] = Element(
    name='Sodium', number=11, symbol='Na', weight=22.99
)
Mg: Final[Element] = Element(
    name='Magnesium', number=12, symbol='Mg', weight=24.305
)
Al: Final[Element] = Element(
    name='Aluminium', number=13, symbol='Al', weight=26.982
)
Si: Final[Element] = Element(
    name='Silicon', number=14, symbol='Si', weight=28.086
)
P: Final[Element] = Element(
    name='Phosphorus', number=15, symbol='P', weight=30.974
)
S: Final[Element] = Element(name='Sulfur', number=16, symbol='S', weight=32.067)
Cl: Final[Element] = Element(
    name='Chlorine', number=17, symbol='Cl', weight=35.453
)
Ar: Final[Element] = Element(
    name='Argon', number=18, symbol='Ar', weight=39.948
)
K: Final[Element] = Element(
    name='Potassium', number=19, symbol='K', weight=39.098
)
Ca: Final[Element] = Element(
    name='Calcium', number=20, symbol='Ca', weight=40.078
)
Sc: Final[Element] = Element(
    name='Scandium', number=21, symbol='Sc', weight=44.956
)
Ti: Final[Element] = Element(
    name='Titanium', number=22, symbol='Ti', weight=47.867
)
V: Final[Element] = Element(
    name='Vanadium', number=23, symbol='V', weight=50.942
)
Cr: Final[Element] = Element(
    name='Chromium', number=24, symbol='Cr', weight=51.996
)
Mn: Final[Element] = Element(
    name='Manganese', number=25, symbol='Mn', weight=54.938
)
Fe: Final[Element] = Element(name='Iron', number=26, symbol='Fe', weight=55.845)
Co: Final[Element] = Element(
    name='Cobalt', number=27, symbol='Co', weight=58.933
)
Ni: Final[Element] = Element(
    name='Nickel', number=28, symbol='Ni', weight=58.693
)
Cu: Final[Element] = Element(
    name='Copper', number=29, symbol='Cu', weight=63.546
)
Zn: Final[Element] = Element(name='Zinc', number=30, symbol='Zn', weight=65.39)
Ga: Final[Element] = Element(
    name='Gallium', number=31, symbol='Ga', weight=69.723
)
Ge: Final[Element] = Element(
    name='Germanium', number=32, symbol='Ge', weight=72.61
)
As: Final[Element] = Element(
    name='Arsenic', number=33, symbol='As', weight=74.922
)
Se: Final[Element] = Element(
    name='Selenium', number=34, symbol='Se', weight=78.96
)
Br: Final[Element] = Element(
    name='Bromine', number=35, symbol='Br', weight=79.904
)
Kr: Final[Element] = Element(
    name='Krypton', number=36, symbol='Kr', weight=83.8
)
Rb: Final[Element] = Element(
    name='Rubidium', number=37, symbol='Rb', weight=85.468
)
Sr: Final[Element] = Element(
    name='Strontium', number=38, symbol='Sr', weight=87.62
)
Y: Final[Element] = Element(
    name='Yttrium', number=39, symbol='Y', weight=88.906
)
Zr: Final[Element] = Element(
    name='Zirconium', number=40, symbol='Zr', weight=91.224
)
Nb: Final[Element] = Element(
    name='Niobiu', number=41, symbol='Nb', weight=92.906
)
Mo: Final[Element] = Element(
    name='Molybdenum', number=42, symbol='Mo', weight=95.94
)
Tc: Final[Element] = Element(
    name='Technetium', number=43, symbol='Tc', weight=98
)
Ru: Final[Element] = Element(
    name='Ruthenium', number=44, symbol='Ru', weight=101.07
)
Rh: Final[Element] = Element(
    name='Rhodium', number=45, symbol='Rh', weight=102.906
)
Pd: Final[Element] = Element(
    name='Palladium', number=46, symbol='Pd', weight=106.42
)
Ag: Final[Element] = Element(
    name='Silver', number=47, symbol='Ag', weight=107.868
)
Cd: Final[Element] = Element(
    name='Cadmium', number=48, symbol='Cd', weight=112.412
)
In: Final[Element] = Element(
    name='Indium', number=49, symbol='In', weight=114.818
)
Sn: Final[Element] = Element(name='Tin', number=50, symbol='Sn', weight=118.711)
Sb: Final[Element] = Element(
    name='Antimony', number=51, symbol='Sb', weight=121.76
)
Te: Final[Element] = Element(
    name='Tellurium', number=52, symbol='Te', weight=127.6
)
I: Final[Element] = Element(
    name='Iodine', number=53, symbol='I', weight=126.904
)
Xe: Final[Element] = Element(
    name='Xenon', number=54, symbol='Xe', weight=131.29
)
Cs: Final[Element] = Element(
    name='Caesium', number=55, symbol='Cs', weight=132.905
)
Ba: Final[Element] = Element(
    name='Barium', number=56, symbol='Ba', weight=137.328
)
La: Final[Element] = Element(
    name='Lanthanum', number=57, symbol='La', weight=138.906
)
Ce: Final[Element] = Element(
    name='Cerium', number=58, symbol='Ce', weight=140.116
)
Pr: Final[Element] = Element(
    name='Praseodymium', number=59, symbol='Pr', weight=140.908
)
Nd: Final[Element] = Element(
    name='Neodymium', number=60, symbol='Nd', weight=144.24
)
Pm: Final[Element] = Element(
    name='Promethium', number=61, symbol='Pm', weight=145
)
Sm: Final[Element] = Element(
    name='Samarium', number=62, symbol='Sm', weight=150.36
)
Eu: Final[Element] = Element(
    name='Europium', number=63, symbol='Eu', weight=151.964
)
Gd: Final[Element] = Element(
    name='Gadolinium', number=64, symbol='Gd', weight=157.25
)
Tb: Final[Element] = Element(
    name='Terbium', number=65, symbol='Tb', weight=158.925
)
Dy: Final[Element] = Element(
    name='Dysprosium', number=66, symbol='Dy', weight=162.5
)
Ho: Final[Element] = Element(
    name='Holmium', number=67, symbol='Ho', weight=164.93
)
Er: Final[Element] = Element(
    name='Erbium', number=68, symbol='Er', weight=167.26
)
Tm: Final[Element] = Element(
    name='Thulium', number=69, symbol='Tm', weight=168.934
)
Yb: Final[Element] = Element(
    name='Ytterbium', number=70, symbol='Yb', weight=173.04
)
Lu: Final[Element] = Element(
    name='Lutetium', number=71, symbol='Lu', weight=174.967
)
Hf: Final[Element] = Element(
    name='Hafnium', number=72, symbol='Hf', weight=178.49
)
Ta: Final[Element] = Element(
    name='Tantalum', number=73, symbol='Ta', weight=180.948
)
W: Final[Element] = Element(
    name='Tungsten', number=74, symbol='W', weight=183.84
)
Re: Final[Element] = Element(
    name='Rhenium', number=75, symbol='Re', weight=186.207
)
Os: Final[Element] = Element(
    name='Osmium', number=76, symbol='Os', weight=190.23
)
Ir: Final[Element] = Element(
    name='Iridium', number=77, symbol='Ir', weight=192.217
)
Pt: Final[Element] = Element(
    name='Platinum', number=78, symbol='Pt', weight=195.078
)
Au: Final[Element] = Element(
    name='Gold', number=79, symbol='Au', weight=196.967
)
Hg: Final[Element] = Element(
    name='Mercury', number=80, symbol='Hg', weight=200.59
)
Tl: Final[Element] = Element(
    name='Thallium', number=81, symbol='Tl', weight=204.383
)
Pb: Final[Element] = Element(name='Lead', number=82, symbol='Pb', weight=207.2)
Bi: Final[Element] = Element(
    name='Bismuth', number=83, symbol='Bi', weight=208.98
)
Po: Final[Element] = Element(
    name='Polonium', number=84, symbol='Po', weight=209
)
At: Final[Element] = Element(
    name='Astatine', number=85, symbol='At', weight=210
)
Rn: Final[Element] = Element(name='Radon', number=86, symbol='Rn', weight=222)
Fr: Final[Element] = Element(
    name='Francium', number=87, symbol='Fr', weight=223
)
Ra: Final[Element] = Element(name='Radium', number=88, symbol='Ra', weight=226)
Ac: Final[Element] = Element(
    name='Actinium', number=89, symbol='Ac', weight=227
)
Th: Final[Element] = Element(
    name='Thorium', number=90, symbol='Th', weight=232.038
)
Pa: Final[Element] = Element(
    name='Protactinium', number=91, symbol='Pa', weight=231.036
)
U: Final[Element] = Element(
    name='Uranium', number=92, symbol='U', weight=238.029
)
Np: Final[Element] = Element(
    name='Neptunium', number=93, symbol='Np', weight=237
)
Pu: Final[Element] = Element(
    name='Plutonium', number=94, symbol='Pu', weight=244
)
Am: Final[Element] = Element(
    name='Americium', number=95, symbol='Am', weight=243
)
Cm: Final[Element] = Element(name='Curium', number=96, symbol='Cm', weight=247)
Bk: Final[Element] = Element(
    name='Berkelium', number=97, symbol='Bk', weight=247
)
Cf: Final[Element] = Element(
    name='Californium', number=98, symbol='Cf', weight=251
)
Es: Final[Element] = Element(
    name='Einsteinium', number=99, symbol='Es', weight=252
)
Fm: Final[Element] = Element(
    name='Fermium', number=100, symbol='Fm', weight=257
)
Md: Final[Element] = Element(
    name='Mendelevium', number=101, symbol='Md', weight=258
)
No: Final[Element] = Element(
    name='Nobelium', number=102, symbol='No', weight=259
)
Lr: Final[Element] = Element(
    name='Lawrencium', number=103, symbol='Lr', weight=262
)
Rf: Final[Element] = Element(
    name='Rutherfordium', number=104, symbol='Rf', weight=267
)
Db: Final[Element] = Element(
    name='Dubnium', number=105, symbol='Db', weight=268
)
Sg: Final[Element] = Element(
    name='Seaborgium', number=106, symbol='Sg', weight=269
)
Bh: Final[Element] = Element(
    name='Bohrium', number=107, symbol='Bh', weight=270
)
Hs: Final[Element] = Element(
    name='Hassium', number=108, symbol='Hs', weight=269
)
Mt: Final[Element] = Element(
    name='Meitnerium', number=109, symbol='Mt', weight=278
)
Ds: Final[Element] = Element(
    name='Darmstadtium', number=110, symbol='Ds', weight=281
)
Rg: Final[Element] = Element(
    name='Roentgenium', number=111, symbol='Rg', weight=281
)
Cn: Final[Element] = Element(
    name='Copernicium', number=112, symbol='Cn', weight=285
)
Nh: Final[Element] = Element(
    name='Nihonium', number=113, symbol='Nh', weight=284
)
Fl: Final[Element] = Element(
    name='Flerovium', number=114, symbol='Fl', weight=289
)
Mc: Final[Element] = Element(
    name='Moscovium', number=115, symbol='Mc', weight=288
)
Lv: Final[Element] = Element(
    name='Livermorium', number=116, symbol='Lv', weight=293
)
Ts: Final[Element] = Element(
    name='Tennessine', number=117, symbol='Ts', weight=292
)
Og: Final[Element] = Element(
    name='Oganesson', number=118, symbol='Og', weight=294
)
# pylint: enable=invalid-name

# fmt: off
# Lanthanides
_L: Final[Sequence[Element]] = (
    La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu)
# Actinides
_A: Final[Sequence[Element]] = (
    Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr)

# pylint: disable=bad-whitespace
PERIODIC_TABLE: Final[Sequence[Element]] = (
    X,  # Unknown
    H,                                                                   He,
    Li, Be,                                          B,  C,  N,  O,  F,  Ne,
    Na, Mg,                                          Al, Si, P,  S,  Cl, Ar,
    K,  Ca,  Sc, Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
    Rb, Sr,  Y,  Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I,  Xe,
    Cs, Ba, *_L, Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, *_A, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og
)
# pylint: enable=bad-whitespace
# fmt: on
ATOMIC_SYMBOL: Mapping[int, str] = {e.number: e.symbol for e in PERIODIC_TABLE}
ATOMIC_NUMBER = {e.symbol: e.number for e in PERIODIC_TABLE}
# Add Deuterium as previous table contained it.
ATOMIC_NUMBER['D'] = 1

ATOMIC_NUMBER: Mapping[str, int] = ATOMIC_NUMBER
ATOMIC_WEIGHT: np.ndarray = np.zeros(len(PERIODIC_TABLE), dtype=np.float64)

for e in PERIODIC_TABLE:
  ATOMIC_WEIGHT[e.number] = e.weight
ATOMIC_WEIGHT.setflags(write=False)
