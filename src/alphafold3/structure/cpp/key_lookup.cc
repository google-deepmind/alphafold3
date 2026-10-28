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
#include "absl/container/flat_hash_set.h"
#include "absl/log/check.h"
#include "absl/types/span.h"

namespace alphafold3 {
namespace {

constexpr size_t kBinarySearchTargetRatio = 32;
constexpr uint64_t kDenseLutRatio = 8;
constexpr uint64_t kByteLutMaxBytes = 64 * 1024;
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

template <typename LookupFn, typename SinkFn>
bool ResolveWithLastKeyCache(absl::Span<const int64_t> query_keys,
                             LookupFn lookup, SinkFn emit) {
  int64_t prev_key = query_keys[0];
  int64_t prev_row = lookup(prev_key);
  if (prev_row < 0) {
    return false;
  }
  emit(0, prev_row);
  for (size_t i = 1; i < query_keys.size(); ++i) {
    if (const int64_t key = query_keys[i]; key != prev_key) {
      prev_row = lookup(key);
      if (prev_row < 0) {
        return false;
      }
      prev_key = key;
    }
    emit(i, prev_row);
  }
  return true;
}

template <typename SinkFn>
bool ResolveKeysToSink(absl::Span<const int64_t> query_keys,
                       absl::Span<const int64_t> table_keys, SinkFn emit) {
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
      emit(i, static_cast<int64_t>(offset));
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
        emit);
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
      emit(i, row);
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
      emit);
}

template <typename IntT, typename PredFn>
void WriteMask(absl::Span<const IntT> query_keys, PredFn pred, bool invert,
               absl::Span<bool> out) {
  for (size_t i = 0; i < query_keys.size(); ++i) {
    out[i] = invert != pred(query_keys[i]);
  }
}

template <typename IntT, typename ContainsFn>
void ContainsWithLastKeyCache(absl::Span<const IntT> query_keys,
                              ContainsFn contains, bool invert,
                              absl::Span<bool> out) {
  IntT prev_key = query_keys[0];
  bool prev_out = invert != contains(prev_key);
  out[0] = prev_out;
  for (size_t i = 1; i < query_keys.size(); ++i) {
    if (const IntT key = query_keys[i]; key != prev_key) {
      prev_out = invert != contains(key);
      prev_key = key;
    }
    out[i] = prev_out;
  }
}

}  // namespace

bool ResolveKeys(absl::Span<const int64_t> query_keys,
                 absl::Span<const int64_t> table_keys,
                 absl::Span<int64_t> out) {
  CHECK_EQ(out.size(), query_keys.size());
  return ResolveKeysToSink(query_keys, table_keys,
                           [out](size_t i, int64_t row) { out[i] = row; });
}

template <typename T>
bool GatherByKey(absl::Span<const int64_t> query_keys,
                 absl::Span<const int64_t> table_keys,
                 absl::Span<const T> table_values, absl::Span<T> out) {
  CHECK_EQ(table_values.size(), table_keys.size());
  CHECK_EQ(out.size(), query_keys.size());
  return ResolveKeysToSink(query_keys, table_keys,
                           [table_values, out](size_t i, int64_t row) {
                             out[i] = table_values[row];
                           });
}

template bool GatherByKey<uint8_t>(absl::Span<const int64_t> query_keys,
                                   absl::Span<const int64_t> table_keys,
                                   absl::Span<const uint8_t> table_values,
                                   absl::Span<uint8_t> out);
template bool GatherByKey<uint16_t>(absl::Span<const int64_t> query_keys,
                                    absl::Span<const int64_t> table_keys,
                                    absl::Span<const uint16_t> table_values,
                                    absl::Span<uint16_t> out);
template bool GatherByKey<uint32_t>(absl::Span<const int64_t> query_keys,
                                    absl::Span<const int64_t> table_keys,
                                    absl::Span<const uint32_t> table_values,
                                    absl::Span<uint32_t> out);
template bool GatherByKey<uint64_t>(absl::Span<const int64_t> query_keys,
                                    absl::Span<const int64_t> table_keys,
                                    absl::Span<const uint64_t> table_values,
                                    absl::Span<uint64_t> out);

template <typename IntT>
void ContainsKeys(absl::Span<const IntT> query_keys,
                  absl::Span<const IntT> table_keys, bool invert,
                  absl::Span<bool> out) {
  using UnsignedT = std::make_unsigned_t<IntT>;
  CHECK_EQ(out.size(), query_keys.size());
  if (query_keys.empty()) {
    return;
  }
  if (table_keys.empty()) {
    absl::c_fill(out, invert);
    return;
  }
  const KeyStats<IntT> stats = AnalyzeKeys(table_keys);
  const uint64_t span = stats.span;
  if (stats.strictly_increasing && span == table_keys.size() - 1) {
    const uint64_t num_rows = table_keys.size();
    WriteMask(
        query_keys,
        [umin = stats.umin, num_rows](IntT key) {
          return uint64_t{static_cast<UnsignedT>(key) - umin} < num_rows;
        },
        invert, out);
    return;
  }
  if (stats.non_decreasing &&
      query_keys.size() <= table_keys.size() / kBinarySearchTargetRatio) {
    ContainsWithLastKeyCache(
        query_keys,
        [table_keys](IntT key) {
          return absl::c_binary_search(table_keys, key);
        },
        invert, out);
    return;
  }
  const uint64_t budget = query_keys.size() + table_keys.size();
  if (span < kByteLutMaxBytes && span / kDenseLutRatio <= budget) {
    std::vector<uint8_t> present(static_cast<size_t>(span) + 1, 0);
    for (const IntT key : table_keys) {
      present[static_cast<size_t>(static_cast<UnsignedT>(key) - stats.umin)] =
          1;
    }
    WriteMask(
        query_keys,
        [umin = stats.umin, span = stats.span,
         table = absl::MakeConstSpan(present)](IntT key) {
          const UnsignedT offset = static_cast<UnsignedT>(key) - umin;
          return offset <= span && table[static_cast<size_t>(offset)] != 0;
        },
        invert, out);
    return;
  }
  if (span < kDenseLutMaxBytes * 8 && (span / 8) / kDenseLutRatio <= budget) {
    std::vector<uint64_t> present(static_cast<size_t>(span >> 6) + 1, 0);
    for (const IntT key : table_keys) {
      const UnsignedT offset = static_cast<UnsignedT>(key) - stats.umin;
      present[static_cast<size_t>(offset >> 6)] |= uint64_t{1} << (offset & 63);
    }
    WriteMask(
        query_keys,
        [umin = stats.umin, span = stats.span,
         bits = absl::MakeConstSpan(present)](IntT key) {
          const UnsignedT offset = static_cast<UnsignedT>(key) - umin;
          return offset <= span &&
                 ((bits[static_cast<size_t>(offset >> 6)] >> (offset & 63)) &
                  1) != 0;
        },
        invert, out);
    return;
  }
  const absl::flat_hash_set<IntT> key_set(table_keys.begin(), table_keys.end());
  ContainsWithLastKeyCache(
      query_keys, [&key_set](IntT key) { return key_set.contains(key); },
      invert, out);
}

template void ContainsKeys<int32_t>(absl::Span<const int32_t> query_keys,
                                    absl::Span<const int32_t> table_keys,
                                    bool invert, absl::Span<bool> out);
template void ContainsKeys<int64_t>(absl::Span<const int64_t> query_keys,
                                    absl::Span<const int64_t> table_keys,
                                    bool invert, absl::Span<bool> out);

}  // namespace alphafold3
