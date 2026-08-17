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

"""Constants shared across modules in the AlphaFold data pipeline."""

from alphafold3.constants import residue_names

MSA_GAP_IDX = residue_names.PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN_AND_GAP.index(
    '-'
)

# Feature groups.
NUM_SEQ_NUM_RES_MSA_FEATURES = ('msa', 'msa_mask', 'deletion_matrix')
NUM_SEQ_MSA_FEATURES = ('msa_species_identifiers',)
TEMPLATE_FEATURES = (
    'template_aatype',
    'template_atom_positions',
    'template_atom_mask',
)
MSA_PAD_VALUES = {'msa': MSA_GAP_IDX, 'msa_mask': 1, 'deletion_matrix': 0}
