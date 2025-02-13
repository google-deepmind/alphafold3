# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

from alphafold3.cpp import cif_dict

class MmcifLayout:
  def atom_range(self, residue_index: int) -> tuple[int, int]: ...
  def chain_starts(self) -> list[int]: ...
  def chains(self) -> list[int]: ...
  def model_offset(self) -> int: ...
  def num_atoms(self) -> int: ...
  def num_chains(self) -> int: ...
  def num_models(self) -> int: ...
  def num_residues(self) -> int: ...
  def residue_range(self, chain_index: int) -> tuple[int, int]: ...
  def residue_starts(self) -> list[int]: ...
  def residues(self) -> list[int]: ...

def from_mmcif(mmcif: cif_dict.CifDict, model_id: str = ...) -> MmcifLayout: ...
