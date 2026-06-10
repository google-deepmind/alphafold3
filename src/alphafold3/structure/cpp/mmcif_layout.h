// Copyright 2024 DeepMind Technologies Limited
//
// AlphaFold 3 source code is licensed under the Apache License, Version 2.0
// (the "License"); you may not use this file except in compliance with the
// License. You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// To request access to the AlphaFold 3 model parameters, follow the process set
// out at https://github.com/google-deepmind/alphafold3. You may only use these
// if received directly from Google. Use is subject to terms of use available at
// https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

#ifndef ALPHAFOLD3_SRC_ALPHAFOLD3_STRUCTURE_PYTHON_MMCIF_LAYOUT_H_
#define ALPHAFOLD3_SRC_ALPHAFOLD3_STRUCTURE_PYTHON_MMCIF_LAYOUT_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "alphafold3/parsers/cpp/cif_dict_lib.h"

namespace alphafold3 {

// Holds the layout of a parsed mmCIF file.
class MmcifLayout {
 public:
  MmcifLayout(std::vector<std::size_t> chain_ends,
              std::vector<std::size_t> residues, std::size_t model_offset,
              std::size_t num_models)
      : chain_ends_(std::move(chain_ends)),
        residue_ends_(std::move(residues)),
        model_offset_(model_offset),
        num_models_(num_models) {}

  // Reads a layout from a valid parsed mmCIF. If a valid model_id is provided
  // the offsets will select that model from the mmCIF.
  // If no model_id is specified, we calculate the layout of the first model
  // only. Therefore it is a requirement that each model has identical atom
  // layouts. An error is returned if the atom counts do not between models.
  static absl::StatusOr<MmcifLayout> Create(const CifDict& mmcif,
                                            absl::string_view model_id = "");

  std::string ToDebugString() const;

  // Returns the start index and one past the last residue index of a given
  // chain. A chain_index of n refers to the n-th chain in the mmCIF. The
  // returned residue indices are 0-based enumerations of residues in the
  // _atom_site records, and therefore do not include missing residues.
  std::pair<std::size_t, std::size_t> residue_range(
      std::size_t chain_index) const {
    if (chain_index > 0) {
      return {chain_ends_[chain_index - 1], chain_ends_[chain_index]};
    } else {
      return {0, chain_ends_[0]};
    }
  }

  // Returns the start index and one past the last index of a given residue.
  // A residue_index of n refers to the n-th residue in the mmCIF, not
  // including residues that are unresolved (i.e. only using _atom_site).
  std::pair<std::size_t, std::size_t> atom_range(
      std::size_t residue_index) const {
    if (residue_index > 0) {
      return {residue_ends_[residue_index - 1], residue_ends_[residue_index]};
    } else {
      return {model_offset_, residue_ends_[residue_index]};
    }
  }

  // If model_id was provided during construction then this is 1, otherwise
  // it is the number of models present in the mmCIF.
  std::size_t num_models() const { return num_models_; }
  // The number of atoms in the chosen model.
  std::size_t num_atoms() const {
    return residue_ends_.empty() ? 0 : residue_ends_.back() - model_offset_;
  }
  // The number of chains in the chosen model.
  std::size_t num_chains() const { return chain_ends_.size(); }
  // The number of residues in the chosen model, not counting unresolved
  // residues.
  std::size_t num_residues() const { return residue_ends_.size(); }

  // Returns the first atom index that is part of the specified chain.
  // The chain is specified using chain_index, which is a 0-based
  // enumeration of the chains in the _atom_site table.
  std::size_t atom_site_from_chain_index(std::size_t chain_index) const {
    if (chain_index == 0) {
      return model_offset_;
    }
    return atom_site_from_residue_index(chain_ends_[chain_index - 1]);
  }

  // Returns the first atom index that is part of the specified residue.
  // The residue is specified using residue_index, which is a 0-based
  // enumeration of the residues in the _atom_site table.
  std::size_t atom_site_from_residue_index(std::size_t residues_index) const {
    if (residues_index == 0) {
      return model_offset_;
    }
    return residue_ends_[residues_index - 1];
  }

  // One past last residue index of each chain. The residue index does not
  // include unresolved residues and is a simple 0-based enumeration of the
  // residues in _atom_site table.
  const std::vector<std::size_t>& chains() const { return chain_ends_; }

  // Indices of the first atom of each chain. Note that this returns atom
  // indices (like residue_starts()), not residue indices (like chains()).
  std::vector<std::size_t> chain_starts() const;

  // One past last atom index of each residue.
  const std::vector<std::size_t>& residues() const { return residue_ends_; }

  // Indices of the first atom of each residue.
  std::vector<std::size_t> residue_starts() const {
    std::vector<std::size_t> residue_starts;
    if (!residue_ends_.empty()) {
      residue_starts.reserve(residue_ends_.size());
      residue_starts.push_back(model_offset_);
      residue_starts.insert(residue_starts.end(), residue_ends_.begin(),
                            residue_ends_.end() - 1);
    }
    return residue_starts;
  }

  // The first atom index that is part of the specified model.
  std::size_t model_offset() const { return model_offset_; }

  void Filter(absl::Span<const std::uint64_t> keep_indices);

 private:
  std::vector<std::size_t> chain_ends_;
  std::vector<std::size_t> residue_ends_;
  std::size_t model_offset_;
  std::size_t num_models_;
};

}  // namespace alphafold3

#endif  // ALPHAFOLD3_SRC_ALPHAFOLD3_STRUCTURE_PYTHON_MMCIF_LAYOUT_H_
