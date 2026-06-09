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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "absl/algorithm/container.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "alphafold3/parsers/cpp/cif_dict_lib.h"
#include "alphafold3/structure/cpp/mmcif_layout.h"

namespace alphafold3 {

std::string MmcifLayout::ToDebugString() const {
  return absl::StrFormat(
      "MmcifLayout(models=%d, chains=%d, num_residues=%d, atoms=%d)",
      num_models(), num_chains(), num_residues(), num_atoms());
}

// Changes layout to match keep_indices removing empty chains/residues.
void MmcifLayout::Filter(absl::Span<const std::uint64_t> keep_indices) {
  if (num_chains() == 0) {
    return;
  }
  // Update residue indices.
  auto keep_it = absl::c_lower_bound(keep_indices, residue_ends_.front());
  for (auto& residue : residue_ends_) {
    while (keep_it != keep_indices.end() && *keep_it < residue) {
      ++keep_it;
    }
    residue = std::distance(keep_indices.begin(), keep_it);
  }
  // Unique residue_ends_ with updating chains.
  auto first = residue_ends_.begin();
  auto tail = first;
  std::size_t num_skipped = 0;
  std::size_t current = 0;
  for (std::size_t& chain_end : chain_ends_) {
    for (auto e = residue_ends_.begin() + chain_end; first != e; ++first) {
      std::size_t next = *first;
      *tail = next;
      if (current != next) {
        current = next;
        ++tail;
      } else {
        ++num_skipped;
      }
    }
    chain_end -= num_skipped;
  }
  residue_ends_.erase(tail, residue_ends_.end());

  current = 0;
  chain_ends_.erase(std::remove_if(chain_ends_.begin(), chain_ends_.end(),
                                   [&current](std::size_t next) {
                                     bool result = current == next;
                                     current = next;
                                     return result;
                                   }),
                    chain_ends_.end());
  model_offset_ = 0;
}

absl::StatusOr<MmcifLayout> MmcifLayout::Create(const CifDict& mmcif,
                                                absl::string_view model_id) {
  auto model_ids = mmcif["_atom_site.pdbx_PDB_model_num"];
  auto chain_ids = mmcif["_atom_site.label_asym_id"];     // chain ID.
  auto label_seq_ids = mmcif["_atom_site.label_seq_id"];  // residue ID.
  auto auth_seq_ids = mmcif["_atom_site.auth_seq_id"];    // author residue ID.
  auto insertion_codes = mmcif["_atom_site.pdbx_PDB_ins_code"];

  if (model_ids.size() != chain_ids.size() ||
      model_ids.size() != label_seq_ids.size() ||
      (model_ids.size() != auth_seq_ids.size() && !auth_seq_ids.empty()) ||
      (model_ids.size() != insertion_codes.size() &&
       !insertion_codes.empty())) {
    return absl::InvalidArgumentError(absl::StrCat(
        "Invalid _atom_site table.",  //
        " len(_atom_site.pdbx_PDB_model_num): ", model_ids.size(),
        " len(_atom_site.label_asym_id): ", chain_ids.size(),
        " len(_atom_site.label_seq_id): ", label_seq_ids.size(),
        " len(_atom_site.auth_seq_id): ", auth_seq_ids.size(),
        " len(_atom_site.pdbx_PDB_ins_code): ", insertion_codes.size()));
  }
  std::size_t num_atoms = model_ids.size();
  if (num_atoms == 0) {
    return MmcifLayout({}, {}, 0, 0);
  }
  std::size_t model_offset = 0;
  std::size_t num_models;
  std::size_t num_atoms_per_model;
  if (model_id.empty()) {
    absl::string_view first_model_id = model_ids.front();

    // Binary search for where the first model ends.
    num_atoms_per_model = std::distance(
        model_ids.begin(),
        absl::c_upper_bound(model_ids, first_model_id, std::not_equal_to<>{}));
    if (num_atoms % num_atoms_per_model != 0) {
      return absl::InvalidArgumentError(absl::StrCat(
          "Each model must have the same number of atoms: (", num_atoms, " % ",
          num_atoms_per_model, " == ", num_atoms % num_atoms_per_model, ")."));
    }
    num_models = num_atoms / num_atoms_per_model;
    // Test boundary conditions for each model hold.
    for (std::size_t i = 1; i < num_models; ++i) {
      if ((model_ids[i * num_atoms_per_model] !=
           model_ids[(i + 1) * num_atoms_per_model - 1]) ||
          (model_ids[i * num_atoms_per_model - 1] ==
           model_ids[i * num_atoms_per_model])) {
        return absl::InvalidArgumentError(
            absl::StrCat("Each model must have the same number of atoms: (",
                         num_atoms, " % ", num_atoms_per_model,
                         " == ", num_atoms % num_atoms_per_model, ")."));
      }
    }
  } else {
    num_models = 1;
    model_offset =
        std::distance(model_ids.begin(), absl::c_find(model_ids, model_id));
    if (model_offset == model_ids.size()) {
      return absl::InvalidArgumentError(
          absl::StrCat("Unknown model_id: ", model_id));
    }
    model_ids.remove_prefix(model_offset);
    chain_ids.remove_prefix(model_offset);
    label_seq_ids.remove_prefix(model_offset);
    if (!auth_seq_ids.empty()) auth_seq_ids.remove_prefix(model_offset);
    if (!insertion_codes.empty()) insertion_codes.remove_prefix(model_offset);

    num_atoms_per_model = std::distance(
        model_ids.begin(), std::upper_bound(model_ids.begin(), model_ids.end(),
                                            model_id, std::not_equal_to<>{}));
    num_atoms = num_atoms_per_model;
  }
  std::vector<std::size_t> residues;
  std::vector<std::size_t> chains;
  absl::string_view chain_id = chain_ids.front();
  if (!auth_seq_ids.empty() && !insertion_codes.empty()) {
    // If author residue IDs are present then these are preferred to
    // label residue IDs because they work for multi-residue ligands (which
    // are given constant "." label residue IDs).
    // NB: Author residue IDs require both the auth_seq_id and the insertion
    // code to be unique.
    absl::string_view auth_seq_id = auth_seq_ids.front();
    absl::string_view insertion_code = insertion_codes.front();
    for (std::size_t i = 1; i < num_atoms_per_model; ++i) {
      if (absl::string_view current_chain_id = chain_ids[i];
          current_chain_id != chain_id) {
        residues.push_back(i + model_offset);
        chains.push_back(residues.size());
        chain_id = current_chain_id;
        auth_seq_id = auth_seq_ids[i];
        insertion_code = insertion_codes[i];
      } else if (absl::string_view current_seq_id = auth_seq_ids[i],
                 current_insertion_code = insertion_codes[i];
                 insertion_code != current_insertion_code ||
                 auth_seq_id != current_seq_id) {
        residues.push_back(i + model_offset);
        auth_seq_id = current_seq_id;
        insertion_code = current_insertion_code;
      }
    }
  } else {
    absl::string_view label_seq_id = label_seq_ids.front();
    for (std::size_t i = 1; i < num_atoms_per_model; ++i) {
      if (absl::string_view current_chain_id = chain_ids[i];
          current_chain_id != chain_id) {
        residues.push_back(i + model_offset);
        chains.push_back(residues.size());
        chain_id = current_chain_id;
        label_seq_id = label_seq_ids[i];
      } else if (absl::string_view current_seq_id = label_seq_ids[i];
                 label_seq_id != current_seq_id) {
        residues.push_back(i + model_offset);
        label_seq_id = current_seq_id;
      }
    }
  }
  residues.push_back(num_atoms_per_model + model_offset);
  chains.push_back(residues.size());
  return MmcifLayout(std::move(chains), std::move(residues), model_offset,
                     num_models);
}

std::vector<std::size_t> MmcifLayout::chain_starts() const {
  std::vector<std::size_t> chain_starts;
  chain_starts.reserve(chain_ends_.size());
  for (std::size_t index = 0; index < chain_ends_.size(); ++index) {
    chain_starts.push_back(atom_site_from_chain_index(index));
  }
  return chain_starts;
}

}  // namespace alphafold3
