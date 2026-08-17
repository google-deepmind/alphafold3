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

"""Names of things in mmCIF format.

See https://www.iucr.org/__data/iucr/cifdic_html/2/cif_mm.dic/index.html
"""

from collections.abc import Mapping, Sequence, Set
from typing import Final

from alphafold3.constants import atom_types
from alphafold3.constants import residue_names


# The following are all possible values for the "_entity.type".
# https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_entity.type.html
BRANCHED_CHAIN: Final[str] = 'branched'
MACROLIDE_CHAIN: Final[str] = 'macrolide'
NON_POLYMER_CHAIN: Final[str] = 'non-polymer'
POLYMER_CHAIN: Final[str] = 'polymer'
WATER: Final[str] = 'water'

CYCLIC_PSEUDO_PEPTIDE_CHAIN: Final[str] = 'cyclic-pseudo-peptide'
DNA_CHAIN: Final[str] = 'polydeoxyribonucleotide'
DNA_RNA_HYBRID_CHAIN: Final[str] = (
    'polydeoxyribonucleotide/polyribonucleotide hybrid'
)
OTHER_CHAIN: Final[str] = 'other'
PEPTIDE_NUCLEIC_ACID_CHAIN: Final[str] = 'peptide nucleic acid'
POLYPEPTIDE_D_CHAIN: Final[str] = 'polypeptide(D)'
PROTEIN_CHAIN: Final[str] = 'polypeptide(L)'
RNA_CHAIN: Final[str] = 'polyribonucleotide'

# Most common _entity_poly.types.
STANDARD_POLYMER_CHAIN_TYPES: Final[Set[str]] = {
    PROTEIN_CHAIN,
    DNA_CHAIN,
    RNA_CHAIN,
}

# Possible values for _entity.type other than polymer and water.
LIGAND_CHAIN_TYPES: Final[Set[str]] = {
    BRANCHED_CHAIN,
    MACROLIDE_CHAIN,
    NON_POLYMER_CHAIN,
}

# Possible values for _entity.type other than polymer.
NON_POLYMER_CHAIN_TYPES: Final[Set[str]] = {
    *LIGAND_CHAIN_TYPES,
    WATER,
}

# Peptide possible values for _entity_poly.type.
PEPTIDE_CHAIN_TYPES: Final[Set[str]] = {
    CYCLIC_PSEUDO_PEPTIDE_CHAIN,
    POLYPEPTIDE_D_CHAIN,
    PROTEIN_CHAIN,
    PEPTIDE_NUCLEIC_ACID_CHAIN,
}


# Nucleic-acid possible values for _entity_poly.type.
NUCLEIC_ACID_CHAIN_TYPES: Final[Set[str]] = {
    RNA_CHAIN,
    DNA_CHAIN,
    DNA_RNA_HYBRID_CHAIN,
}

# All possible values for _entity_poly.type.
POLYMER_CHAIN_TYPES: Final[Set[str]] = {
    *NUCLEIC_ACID_CHAIN_TYPES,
    *PEPTIDE_CHAIN_TYPES,
    OTHER_CHAIN,
}


TERMINAL_OXYGENS: Final[Mapping[str, str]] = {
    PROTEIN_CHAIN: 'OXT',
    DNA_CHAIN: 'OP3',
    RNA_CHAIN: 'OP3',
}


# For each chain type, which atom should be used to represent each residue.
RESIDUE_REPRESENTATIVE_ATOMS: Final[Mapping[str, str]] = {
    PROTEIN_CHAIN: atom_types.CA,
    DNA_CHAIN: atom_types.C1PRIME,
    RNA_CHAIN: atom_types.C1PRIME,
}

# Methods involving crystallization. See the documentation at
# mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_exptl.method.html
# for the full list of experimental methods.
CRYSTALLIZATION_METHODS: Final[Set[str]] = {
    'X-RAY DIFFRACTION',
    'NEUTRON DIFFRACTION',
    'ELECTRON CRYSTALLOGRAPHY',
    'POWDER CRYSTALLOGRAPHY',
    'FIBER DIFFRACTION',
}

# Possible bond types.
COVALENT_BOND: Final[str] = 'covale'
HYDROGEN_BOND: Final[str] = 'hydrog'
METAL_COORDINATION: Final[str] = 'metalc'
DISULFIDE_BRIDGE: Final[str] = 'disulf'


def is_standard_polymer_type(chain_type: str) -> bool:
  """Returns if chain type is a protein, DNA or RNA chain type.

  Args:
     chain_type: The type of the chain.

  Returns:
    A bool for if the chain_type matches protein, DNA, or RNA.
  """
  return chain_type in STANDARD_POLYMER_CHAIN_TYPES


def guess_polymer_type(chain_residues: Sequence[str]) -> str:
  """Guess the polymer type (protein/rna/dna/other) based on the residues.

  The polymer type is guessed by first checking for any of the standard
  protein residues. If one is present then the chain is considered to be a
  polypeptide. Otherwise we decide by counting residue types and deciding by
  majority voting (e.g. mostly DNA residues -> DNA). If there is a tie between
  the counts, the ordering is rna > dna > other.

  Note that we count MSE and UNK as protein residues.

  Args:
    chain_residues: A sequence of full residue name (1-letter for DNA, 2-letters
      for RNA, 3 for protein). The _atom_site.label_comp_id column in mmCIF.

  Returns:
    The most probable chain type as set in the _entity_poly mmCIF table:
    protein - polypeptide(L), rna - polyribonucleotide,
    dna - polydeoxyribonucleotide or other.
  """
  residue_types = {
      **{r: RNA_CHAIN for r in residue_names.RNA_TYPES},
      **{r: DNA_CHAIN for r in residue_names.DNA_TYPES},
      **{r: PROTEIN_CHAIN for r in residue_names.PROTEIN_TYPES_WITH_UNKNOWN},
      residue_names.MSE: PROTEIN_CHAIN,
  }

  counts = {PROTEIN_CHAIN: 0, RNA_CHAIN: 0, DNA_CHAIN: 0, OTHER_CHAIN: 0}
  for residue in chain_residues:
    residue_type = residue_types.get(residue, OTHER_CHAIN)
    # If we ever see a protein residue we'll consider this a polypeptide(L).
    if residue_type == PROTEIN_CHAIN:
      return residue_type
    counts[residue_type] += 1

  # Make sure protein > rna > dna > other if there is a tie.
  tie_braker = {PROTEIN_CHAIN: 3, RNA_CHAIN: 2, DNA_CHAIN: 1, OTHER_CHAIN: 0}

  def order_fn(item):
    name, count = item
    return count, tie_braker[name]

  most_probable_type = max(counts.items(), key=order_fn)[0]
  return most_probable_type


def fix_non_standard_polymer_res(*, res_name: str, chain_type: str) -> str:
  """Returns the res_name of the closest standard protein/RNA/DNA residue.

  Optimized for the case where a single residue needs to be converted.

  If res_name is already a standard type, it is returned unaltered.
  If a match cannot be found, returns 'UNK' for protein chains and 'N' for
    RNA/DNA chains.

  Args:
     res_name: A residue_name (monomer code from the CCD).
     chain_type: The type of the chain, must be PROTEIN_CHAIN, RNA_CHAIN or
       DNA_CHAIN.

  Returns:
    An element from PROTEIN_TYPES_WITH_UNKNOWN | RNA_TYPES | DNA_TYPES | {'N'}.

  Raises:
    ValueError: If chain_type not in PEPTIDE_CHAIN_TYPES or
      {OTHER_CHAIN, RNA_CHAIN, DNA_CHAIN, DNA_RNA_HYBRID_CHAIN}.
  """
  # Map to one letter code, then back to common res_names.
  one_letter_code = residue_names.letters_three_to_one(res_name, default='X')

  if chain_type in PEPTIDE_CHAIN_TYPES or chain_type == OTHER_CHAIN:
    return residue_names.PROTEIN_COMMON_ONE_TO_THREE.get(one_letter_code, 'UNK')
  elif chain_type == RNA_CHAIN:
    # RNA's CCD monomer code is single-letter.
    return (
        one_letter_code if one_letter_code in residue_names.RNA_TYPES else 'N'
    )
  elif chain_type == DNA_CHAIN:
    return residue_names.DNA_COMMON_ONE_TO_TWO.get(one_letter_code, 'N')
  elif chain_type == DNA_RNA_HYBRID_CHAIN:
    return (
        res_name
        if res_name in residue_names.NUCLEIC_TYPES_WITH_UNKNOWN
        else 'N'
    )
  else:
    raise ValueError(f'Expected a protein/DNA/RNA chain but got {chain_type}')
