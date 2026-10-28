// Copyright 2026 DeepMind Technologies Limited
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

#ifndef ALPHAFOLD3_SRC_ALPHAFOLD3_STRUCTURE_CPP_KEY_LOOKUP_H_
#define ALPHAFOLD3_SRC_ALPHAFOLD3_STRUCTURE_CPP_KEY_LOOKUP_H_

#include <cstdint>

#include "absl/types/span.h"

namespace alphafold3 {

// Resolves each element of `query_keys` to its 0-based position in `table_keys`
// (last occurrence wins on duplicates) and writes `out[i]`. Returns false
// without throwing if any key is absent.
bool ResolveKeys(absl::Span<const int64_t> query_keys,
                 absl::Span<const int64_t> table_keys, absl::Span<int64_t> out);

// Gathers `table_values[row]` for each key in `query_keys` into `out`.
template <typename T>
bool GatherByKey(absl::Span<const int64_t> query_keys,
                 absl::Span<const int64_t> table_keys,
                 absl::Span<const T> table_values, absl::Span<T> out);

extern template bool GatherByKey<uint8_t>(
    absl::Span<const int64_t> query_keys, absl::Span<const int64_t> table_keys,
    absl::Span<const uint8_t> table_values, absl::Span<uint8_t> out);
extern template bool GatherByKey<uint16_t>(
    absl::Span<const int64_t> query_keys, absl::Span<const int64_t> table_keys,
    absl::Span<const uint16_t> table_values, absl::Span<uint16_t> out);
extern template bool GatherByKey<uint32_t>(
    absl::Span<const int64_t> query_keys, absl::Span<const int64_t> table_keys,
    absl::Span<const uint32_t> table_values, absl::Span<uint32_t> out);
extern template bool GatherByKey<uint64_t>(
    absl::Span<const int64_t> query_keys, absl::Span<const int64_t> table_keys,
    absl::Span<const uint64_t> table_values, absl::Span<uint64_t> out);

// Writes `invert != (query_keys[i] in table_keys)` into `out[i]`.
template <typename IntT>
void ContainsKeys(absl::Span<const IntT> query_keys,
                  absl::Span<const IntT> table_keys, bool invert,
                  absl::Span<bool> out);

extern template void ContainsKeys<int32_t>(absl::Span<const int32_t> query_keys,
                                           absl::Span<const int32_t> table_keys,
                                           bool invert, absl::Span<bool> out);
extern template void ContainsKeys<int64_t>(absl::Span<const int64_t> query_keys,
                                           absl::Span<const int64_t> table_keys,
                                           bool invert, absl::Span<bool> out);

}  // namespace alphafold3

#endif  // ALPHAFOLD3_SRC_ALPHAFOLD3_STRUCTURE_CPP_KEY_LOOKUP_H_
