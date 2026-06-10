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

#include <cstdint>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "pybind11/cast.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11_abseil/absl_casters.h"

namespace {

namespace py = pybind11;

absl::flat_hash_map<int64_t, std::vector<int64_t>> IndicesGroupedByValue(
    absl::Span<const int64_t> values) {
  absl::flat_hash_map<int64_t, std::vector<int64_t>> group_indices;
  for (int64_t i = 0, e = values.size(); i < e; ++i) {
    group_indices[values[i]].push_back(i);
  }
  return group_indices;
}

constexpr char kIndicesGroupedByValue[] = R"(
Returns a map from value to a list of indices this value occupies.

E.g. indices_grouped_by_value([1, 1, 2, 3, 3, 1, 1]) returns:
{1: [0, 1, 5, 6], 2: [2], 3: [3, 4]}

Args:
  values: a list of values to group.
)";

}  // namespace

namespace alphafold3 {

void RegisterModuleAggregation(py::module m) {
  m.def("indices_grouped_by_value", &IndicesGroupedByValue, py::arg("values"),
        py::doc(kIndicesGroupedByValue + 1),
        py::call_guard<py::gil_scoped_release>());
}

}  // namespace alphafold3
