# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Sets of chemical components."""

from alphafold3.common import resources
from alphafold3.common import safe_pickle
from typing import Final


_CCD_SETS_CCD_PICKLE_FILE = resources.filename(
    resources.ROOT / 'constants/converters/chemical_component_sets.pickle'
)

with open(_CCD_SETS_CCD_PICKLE_FILE, 'rb') as _f:
  _CCD_SET = safe_pickle.safe_load(_f)

# Glycan (or 'Saccharide') ligands.
# _chem_comp.type containing 'saccharide' and 'linking' (when lower-case).
GLYCAN_LINKING_LIGANDS: Final[frozenset[str]] = _CCD_SET['glycans_linking']

# _chem_comp.type containing 'saccharide' and not 'linking' (when lower-case).
GLYCAN_OTHER_LIGANDS: Final[frozenset[str]] = _CCD_SET['glycans_other']

# Each of these molecules appears in over 1k PDB structures, are used to
# facilitate crystallization conditions, but do not have biological relevance.
COMMON_CRYSTALLIZATION_AIDS: Final[frozenset[str]] = frozenset({
    'SO4', 'GOL', 'EDO', 'PO4', 'ACT', 'PEG', 'DMS', 'TRS', 'PGE', 'PG4', 'FMT',
    'EPE', 'MPD', 'MES', 'CD', 'IOD',
})  # pyformat: disable
