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

"""Batch dataclass."""
import dataclasses
from typing import Self

from alphafold3.model import features
import jax


@dataclasses.dataclass(frozen=True)
class Batch:
  """Dataclass containing batch."""

  msa: features.MSA
  templates: features.Templates
  token_features: features.TokenFeatures
  ref_structure: features.RefStructure
  predicted_structure_info: features.PredictedStructureInfo
  polymer_ligand_bond_info: features.PolymerLigandBondInfo
  ligand_ligand_bond_info: features.LigandLigandBondInfo
  pseudo_beta_info: features.PseudoBetaInfo
  atom_cross_att: features.AtomCrossAtt
  convert_model_output: features.ConvertModelOutput
  frames: features.Frames

  @property
  def num_res(self) -> int:
    return self.token_features.aatype.shape[-1]

  @classmethod
  def from_data_dict(cls, batch: features.BatchDict) -> Self:
    """Construct batch object from dictionary."""
    return cls(
        msa=features.MSA.from_data_dict(batch),
        templates=features.Templates.from_data_dict(batch),
        token_features=features.TokenFeatures.from_data_dict(batch),
        ref_structure=features.RefStructure.from_data_dict(batch),
        predicted_structure_info=features.PredictedStructureInfo.from_data_dict(
            batch
        ),
        polymer_ligand_bond_info=features.PolymerLigandBondInfo.from_data_dict(
            batch
        ),
        ligand_ligand_bond_info=features.LigandLigandBondInfo.from_data_dict(
            batch
        ),
        pseudo_beta_info=features.PseudoBetaInfo.from_data_dict(batch),
        atom_cross_att=features.AtomCrossAtt.from_data_dict(batch),
        convert_model_output=features.ConvertModelOutput.from_data_dict(batch),
        frames=features.Frames.from_data_dict(batch),
    )

  def as_data_dict(self) -> features.BatchDict:
    """Converts batch object to dictionary."""
    output = {
        **self.msa.as_data_dict(),
        **self.templates.as_data_dict(),
        **self.token_features.as_data_dict(),
        **self.ref_structure.as_data_dict(),
        **self.predicted_structure_info.as_data_dict(),
        **self.polymer_ligand_bond_info.as_data_dict(),
        **self.ligand_ligand_bond_info.as_data_dict(),
        **self.pseudo_beta_info.as_data_dict(),
        **self.atom_cross_att.as_data_dict(),
        **self.convert_model_output.as_data_dict(),
        **self.frames.as_data_dict(),
    }
    return output


jax.tree_util.register_dataclass(
    Batch,
    data_fields=[f.name for f in dataclasses.fields(Batch)],
    meta_fields=[],
)
