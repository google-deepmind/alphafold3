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

#include "alphafold3/structure/cpp/key_lookup.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <type_traits>
#include <vector>

#include "absl/algorithm/container.h"
#include "absl/container/flat_hash_map.h"
#include "absl/log/check.h"
#include "absl/types/span.h"

namespace alphafold3 {
namespace {

constexpr size_t kBinarySearchTargetRatio = 32;
constexpr uint64_t kDenseLutRatio = 8;
constexpr uint64_t kDenseLutMaxBytes = 16 * 1024 * 1024;

template <typename IntT>
struct KeyStats {
  using UnsignedT = std::make_unsigned_t<IntT>;
  UnsignedT umin = 0;
  UnsignedT span = 0;
  bool non_decreasing = true;
  bool strictly_increasing = true;
};

template <typename IntT>
KeyStats<IntT> AnalyzeKeys(absl::Span<const IntT> table_keys) {
  using UnsignedT = std::make_unsigned_t<IntT>;
  IntT min_key = table_keys[0];
  IntT max_key = min_key;
  bool non_decreasing = true;
  bool strictly_increasing = true;
  for (size_t i = 1; i < table_keys.size(); ++i) {
    const IntT key = table_keys[i];
    non_decreasing &= (key >= table_keys[i - 1]);
    strictly_increasing &= (key > table_keys[i - 1]);
    min_key = std::min(min_key, key);
    max_key = std::max(max_key, key);
  }
  const UnsignedT umin = static_cast<UnsignedT>(min_key);
  return {umin, static_cast<UnsignedT>(max_key) - umin, non_decreasing,
          strictly_increasing};
}

template <typename LookupFn>
bool ResolveWithLastKeyCache(absl::Span<const int64_t> query_keys,
                             LookupFn lookup, absl::Span<int64_t> out) {
  int64_t prev_key = query_keys[0];
  int64_t prev_row = lookup(prev_key);
  if (prev_row < 0) {
    return false;
  }
  out[0] = prev_row;
  for (size_t i = 1; i < query_keys.size(); ++i) {
    if (const int64_t key = query_keys[i]; key != prev_key) {
      prev_row = lookup(key);
      if (prev_row < 0) {
        return false;
      }
      prev_key = key;
    }
    out[i] = prev_row;
  }
  return true;
}

}  // namespace

bool ResolveKeys(absl::Span<const int64_t> query_keys,
                 absl::Span<const int64_t> table_keys,
                 absl::Span<int64_t> out) {
  CHECK_EQ(out.size(), query_keys.size());
  if (query_keys.empty()) {
    return true;
  }
  if (table_keys.empty()) {
    return false;
  }
  const KeyStats<int64_t> stats = AnalyzeKeys(table_keys);
  if (stats.strictly_increasing && stats.span == table_keys.size() - 1) {
    const uint64_t num_rows = table_keys.size();
    for (size_t i = 0; i < query_keys.size(); ++i) {
      const uint64_t offset = static_cast<uint64_t>(query_keys[i]) - stats.umin;
      if (offset >= num_rows) {
        return false;
      }
      out[i] = static_cast<int64_t>(offset);
    }
    return true;
  }
  if (stats.non_decreasing &&
      query_keys.size() <= table_keys.size() / kBinarySearchTargetRatio) {
    return ResolveWithLastKeyCache(
        query_keys,
        [&](int64_t key) -> int64_t {
          auto it = absl::c_upper_bound(table_keys, key);
          return (it != table_keys.begin() && *std::prev(it) == key)
                     ? static_cast<int64_t>(std::prev(it) - table_keys.begin())
                     : -1;
        },
        out);
  }
  if (table_keys.size() <=
          static_cast<size_t>(std::numeric_limits<int32_t>::max()) &&
      stats.span < kDenseLutMaxBytes / sizeof(int32_t) &&
      (stats.span * sizeof(int32_t)) / kDenseLutRatio <=
          query_keys.size() + table_keys.size()) {
    std::vector<int32_t> lut(static_cast<size_t>(stats.span) + 1, -1);
    for (size_t i = 0; i < table_keys.size(); ++i) {
      const uint64_t offset = static_cast<uint64_t>(table_keys[i]) - stats.umin;
      if (offset > stats.span) {
        return false;
      }
      lut[static_cast<size_t>(offset)] = static_cast<int32_t>(i);
    }
    for (size_t i = 0; i < query_keys.size(); ++i) {
      const uint64_t offset = static_cast<uint64_t>(query_keys[i]) - stats.umin;
      if (offset > stats.span) {
        return false;
      }
      const int32_t row = lut[static_cast<size_t>(offset)];
      if (row < 0) {
        return false;
      }
      out[i] = row;
    }
    return true;
  }
  absl::flat_hash_map<int64_t, int64_t> index_by_key;
  index_by_key.reserve(table_keys.size());
  for (size_t i = 0; i < table_keys.size(); ++i) {
    const int64_t key = table_keys[i];
    index_by_key[key] = static_cast<int64_t>(i);
  }
  return ResolveWithLastKeyCache(
      query_keys,
      [&](int64_t key) -> int64_t {
        auto it = index_by_key.find(key);
        return it != index_by_key.end() ? it->second : -1;
      },
      out);
}

}  // namespace alphafold3
