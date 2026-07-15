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
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "pybind11/cast.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11_abseil/absl_casters.h"

namespace {

namespace py = pybind11;

py::array_t<bool> IsIn(const py::array_t<int64_t, py::array::c_style>& array,
                       const absl::flat_hash_set<int64_t>& test_elements,
                       bool invert) {
  const size_t num_elements = array.size();

  py::array_t<bool> output(num_elements);
  std::fill(output.mutable_data(), output.mutable_data() + output.size(),
            invert);

  // Shortcut: The output will be trivially always false if test_elements empty.
  if (test_elements.empty()) {
    return output;
  }

  for (size_t i = 0; i < num_elements; ++i) {
    if (test_elements.contains(array.data()[i])) {
      output.mutable_data()[i] = !invert;
    }
  }
  if (array.ndim() > 1) {
    auto shape =
        std::vector<ptrdiff_t>(array.shape(), array.shape() + array.ndim());
    return output.reshape(shape);
  }
  return output;
}

constexpr char kIsInDoc[] = R"(
Computes whether each element is in test_elements.

Same use as np.isin, but much faster. If len(array) = n, len(test_elements) = m:
* This function has complexity O(n).
* np.isin with kind='sort' has complexity O(m*log(m) + n * log(m)).

Args:
  array: Input NumPy array with dtype=np.int64.
  test_elements: The values against which to test each value of array.
  invert: If True, the values in the returned array are inverted, as if
    calculating `element not in test_elements`. Default is False.
    `isin(a, b, invert=True)` is equivalent to but faster than `~isin(a, b)`.

Returns
  A boolean array of the same shape as the input array. Each value `val` is:
  * `val in test_elements` if `invert=False`,
  * `val not in test_elements` if `invert=True`.
)";

}  // namespace

namespace alphafold3 {

void RegisterModuleMembership(pybind11::module m) {
  m.def("isin", &IsIn, py::arg("array"), py::arg("test_elements"),
        py::kw_only(), py::arg("invert") = false, py::doc(kIsInDoc + 1));
}

}  // namespace alphafold3
