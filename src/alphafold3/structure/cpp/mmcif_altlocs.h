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

#ifndef ALPHAFOLD3_SRC_ALPHAFOLD3_STRUCTURE_PYTHON_MMCIF_ALTLOCS_H_
#define ALPHAFOLD3_SRC_ALPHAFOLD3_STRUCTURE_PYTHON_MMCIF_ALTLOCS_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "absl/types/span.h"
#include "alphafold3/structure/cpp/mmcif_layout.h"

namespace alphafold3 {

// Returns the list of indices that should be kept after resolving alt-locs.
// 1) Partial Residue. Each cycle of alt-locs are resolved separately with the
//    highest occupancy alt-loc. Tie-breaks are resolved alphabetically. See
//    tests for examples.
// 2) Whole Residue. These are resolved in two passes.
//    a) The residue with the highest occupancy is chosen.
//    b) The locations for a given residue are resolved.
//    All tie-breaks are resolved alphabetically. See tests for examples.
//
// Preconditions: layout and comp_ids, alt_ids, occupancies are all from same
// mmCIF file and chain_indices are monotonically increasing and less than
// layout.num_chains().
//
// comp_ids from '_atom_site.label_comp_id'.
// alt_ids from '_atom_site.label_alt_id'.
// occupancies from '_atom_site.occupancy'.
std::vector<std::uint64_t> ResolveMmcifAltLocs(
    const MmcifLayout& layout, absl::Span<const std::string> comp_ids,
    absl::Span<const std::string> atom_ids,
    absl::Span<const std::string> alt_ids,
    absl::Span<const std::string> occupancies,
    absl::Span<const std::size_t> chain_indices);

}  // namespace alphafold3

#endif  // ALPHAFOLD3_SRC_ALPHAFOLD3_STRUCTURE_PYTHON_MMCIF_ALTLOCS_H_
