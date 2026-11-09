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

#include "alphafold3/structure/cpp/mmcif_altlocs.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "absl/algorithm/container.h"
#include "absl/log/log.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "alphafold3/structure/cpp/mmcif_layout.h"

namespace alphafold3 {
namespace {

float OccupancyToFloat(absl::string_view occupancy) {
  float result = 0.0f;
  LOG_IF(ERROR, !absl::SimpleAtof(occupancy, &result))
      << "Invalid Occupancy: " << occupancy;
  return result;
}

// Deuterium is the same atom as Hydrogen so keep equivalent for grouping.
bool AtomEquiv(absl::string_view lhs, absl::string_view rhs) {
  if (lhs == rhs) return true;
  if (lhs.empty() != rhs.empty()) return false;
  // Both lhs and rhs are guaranteed to be non-empty after this.
  char first_lhs = lhs.front();
  char second_rhs = rhs.front();
  if ((first_lhs == 'H' && second_rhs == 'D') ||
      (first_lhs == 'D' && second_rhs == 'H')) {
    lhs.remove_prefix(1);
    rhs.remove_prefix(1);
    return lhs == rhs;
  }
  return false;
}

// Calls group_callback with that start index and count for each group of
// equivalent values in `values`, starting at `start` and ending at `count`.
// Example:
// GroupBy({"B", "B", "B", "C", "C"}, 0, 5, [](size_t start, size_t count) {
//   absl::Printf("start=%d, count=%d\n", start, count);
// });
// Would print:
// start=0, count=3
// start=3, count=2
template <typename GroupCallback,
          typename IsEqual = std::equal_to<absl::string_view>>
void GroupBy(absl::Span<const std::string> values, std::size_t start,
             std::size_t count, GroupCallback&& group_callback,
             IsEqual&& is_equal = std::equal_to<absl::string_view>{}) {
  std::size_t span_start = start;
  if (count > 0) {
    for (std::size_t i = start + 1; i < start + count; ++i) {
      if (!is_equal(values[i], values[span_start])) {
        group_callback(span_start, i - span_start);
        span_start = i;
      }
    }
    group_callback(span_start, start + count - span_start);
  }
}

void ProcessAltLocGroupsWhole(std::size_t alt_loc_start,
                              std::size_t alt_loc_count,
                              absl::Span<const std::string> comp_ids,
                              absl::Span<const std::string> atom_ids,
                              absl::Span<const std::string> alt_ids,
                              absl::Span<const std::string> occupancies,
                              std::vector<std::uint64_t>& in_out_keep_indices) {
  std::pair<std::size_t, std::size_t> best_split = {alt_loc_start,
                                                    alt_loc_count};
  std::vector<absl::string_view> alt_loc_groups;
  float best_occupancy = -std::numeric_limits<float>::infinity();
  absl::string_view best_group = alt_ids[alt_loc_start];
  std::vector<std::pair<std::size_t, float>> occupancy_stats;

  // Group by residue type.
  GroupBy(comp_ids, alt_loc_start, alt_loc_count,
          [&](std::size_t start, std::size_t count) {
            // This callback selects the best residue group and the best
            // alt-loc ID within that group.
            alt_loc_groups.clear();
            occupancy_stats.clear();
            // Calculate total occupancy for residue type.
            for (std::size_t i = 0; i < count; ++i) {
              absl::string_view alt_loc_id = alt_ids[start + i];
              float occupancy = OccupancyToFloat(occupancies[start + i]);
              if (auto loc = absl::c_find(alt_loc_groups, alt_loc_id);
                  loc == alt_loc_groups.end()) {
                occupancy_stats.emplace_back(1, occupancy);
                alt_loc_groups.push_back(alt_loc_id);
              } else {
                auto& stat =
                    occupancy_stats[std::distance(alt_loc_groups.begin(), loc)];
                ++stat.first;
                stat.second += occupancy;
              }
            }
            float total_occupancy = 0.0;
            for (auto& stat : occupancy_stats) {
              total_occupancy += stat.second / stat.first;
            }
            absl::string_view group = *absl::c_min_element(alt_loc_groups);
            // Compares occupancy of residue to best seen so far.
            // Tie breaks alphabetic.
            if (total_occupancy > best_occupancy ||
                (total_occupancy == best_occupancy && group < best_group)) {
              // Selects the best sub group.
              best_group = alt_loc_groups.front();
              float best_amount = occupancy_stats.front().second /
                                  occupancy_stats.front().first;
              for (std::size_t i = 1; i < occupancy_stats.size(); ++i) {
                float amount =
                    occupancy_stats[i].second / occupancy_stats[i].first;
                if (amount > best_amount ||
                    (amount == best_amount && alt_loc_groups[i] < best_group)) {
                  best_amount = amount;
                  best_group = alt_loc_groups[i];
                }
              }
              best_occupancy = total_occupancy;
              best_split = {start, count};
            }
          });

  // Now that the best residue type has been selected and the best alt-loc
  // within that has been selected add indices of indices to keep to the keep
  // list.
  auto [split_start, split_count] = best_split;
  GroupBy(
      atom_ids, split_start, split_count,
      [&in_out_keep_indices, &alt_ids, best_group](std::size_t start,
                                                   std::size_t count) {
        // This makes sure we select an atom for each atom id even if it does
        // not have our selected alt-loc ID.
        std::size_t best_index = start;
        for (std::size_t i = 0; i < count; ++i) {
          if (alt_ids[start + i] == best_group) {
            best_index = start + i;
            break;
          }
        }
        in_out_keep_indices.push_back(best_index);
      },
      AtomEquiv);
}

// Finds the alt-loc group with the highest score and pushes the indices on to
// the back of in_out_keep_indices.
void ProcessAltLocGroupPartial(
    std::size_t alt_loc_start, std::size_t alt_loc_count,
    absl::Span<const std::string> atom_ids,
    absl::Span<const std::string> alt_ids,
    absl::Span<const std::string> occupancies,
    std::vector<std::uint64_t>& in_out_keep_indices) {
  GroupBy(
      atom_ids, alt_loc_start, alt_loc_count,
      [&](std::size_t start, std::size_t count) {
        if (count == 1) {
          in_out_keep_indices.push_back(start);
        } else {
          float best_occ = OccupancyToFloat(occupancies[start]);
          std::size_t best_index = start;
          absl::string_view best_group = alt_ids[start];
          for (std::size_t i = 1; i < count; ++i) {
            float occ = OccupancyToFloat(occupancies[start + i]);
            absl::string_view group = alt_ids[start + i];
            if (occ > best_occ || (occ == best_occ && group < best_group)) {
              best_group = group;
              best_index = start + i;
              best_occ = occ;
            }
          }
          in_out_keep_indices.push_back(best_index);
        }
      },
      AtomEquiv);
}

}  // namespace

// Resolves alt-locs returning the atom indices that will be left.
absl::StatusOr<std::vector<std::uint64_t>> ResolveMmcifAltLocs(
    const MmcifLayout& layout, absl::Span<const std::string> comp_ids,
    absl::Span<const std::string> atom_ids,
    absl::Span<const std::string> alt_ids,
    absl::Span<const std::string> occupancies,
    absl::Span<const std::size_t> chain_indices) {
  const std::size_t required_atoms = layout.model_offset() + layout.num_atoms();
  if (comp_ids.size() < required_atoms || atom_ids.size() < required_atoms ||
      alt_ids.size() < required_atoms || occupancies.size() < required_atoms) {
    return absl::InvalidArgumentError(absl::StrCat(
        "Invalid alt-loc input sizes: required >= ", required_atoms,
        ", got comp_ids=", comp_ids.size(), ", atom_ids=", atom_ids.size(),
        ", alt_ids=", alt_ids.size(), ", occupancies=", occupancies.size()));
  }

  std::vector<std::uint64_t> keep_indices;
  keep_indices.reserve(layout.num_atoms());
  std::size_t alt_loc_start = 0;
  for (std::size_t chain_index : chain_indices) {
    if (chain_index >= layout.num_chains()) {
      return absl::InvalidArgumentError(
          absl::StrCat("Chain index ", chain_index, " out of range for ",
                       layout.num_chains(), " chains."));
    }
    auto [residues_start, residues_end] = layout.residue_range(chain_index);
    for (std::size_t residue = residues_start; residue < residues_end;
         ++residue) {
      std::size_t alt_loc_count = 0;
      auto [atom_start, atom_end] = layout.atom_range(residue);
      for (std::size_t i = atom_start; i < atom_end; ++i) {
        absl::string_view alt_loc_id = alt_ids[i];
        if (alt_loc_id.empty()) {
          return absl::InvalidArgumentError(
              "Invalid '_atom_site.label_alt_id': contains empty string.");
        }
        if (alt_loc_id == "." || alt_loc_id == "?") {
          if (alt_loc_count > 0) {
            ProcessAltLocGroupPartial(alt_loc_start, alt_loc_count, atom_ids,
                                      alt_ids, occupancies, keep_indices);
            alt_loc_count = 0;
          }
          keep_indices.push_back(i);
        } else {
          if (alt_loc_count == 0) {
            alt_loc_start = i;
          }
          ++alt_loc_count;
        }
      }
      if (alt_loc_count > 0) {
        if (atom_end - atom_start == alt_loc_count) {
          ProcessAltLocGroupsWhole(alt_loc_start, alt_loc_count, comp_ids,
                                   atom_ids, alt_ids, occupancies,
                                   keep_indices);
        } else {
          ProcessAltLocGroupPartial(alt_loc_start, alt_loc_count, atom_ids,
                                    alt_ids, occupancies, keep_indices);
        }
      }
    }
  }

  return keep_indices;
}

}  // namespace alphafold3
