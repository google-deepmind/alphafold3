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

#include "absl/strings/str_cat.h"
#include "pybind11/cast.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"

namespace {

namespace py = pybind11;

py::array_t<double> ComputeMsaProfile(
    const py::array_t<int, py::array::c_style>& msa, int num_residue_types) {
  if (msa.size() == 0) {
    throw py::value_error("The MSA must be non-empty.");
  }
  if (msa.ndim() != 2) {
    throw py::value_error(absl::StrCat("The MSA must be rectangular, got ",
                                       msa.ndim(), "-dimensional MSA array."));
  }
  const int msa_depth = msa.shape()[0];
  const int sequence_length = msa.shape()[1];

  py::array_t<double> profile({sequence_length, num_residue_types});
  std::fill(profile.mutable_data(), profile.mutable_data() + profile.size(),
            0.0f);
  auto profile_unchecked = profile.mutable_unchecked<2>();

  const double normalized_count = 1.0 / msa_depth;
  const int* msa_it = msa.data();
  for (int row_index = 0; row_index < msa_depth; ++row_index) {
    for (int column_index = 0; column_index < sequence_length; ++column_index) {
      const int residue_code = *(msa_it++);
      if (residue_code < 0 || residue_code >= num_residue_types) {
        throw py::value_error(
            absl::StrCat("All residue codes must be positive and smaller than "
                         "num_residue_types ",
                         num_residue_types, ", got ", residue_code));
      }
      profile_unchecked(column_index, residue_code) += normalized_count;
    }
  }
  return profile;
}

constexpr char kComputeMsaProfileDoc[] = R"(
Computes MSA profile for the given encoded MSA.

Args:
  msa: A Numpy array of shape (num_msa, num_res) with the integer coded MSA.
  num_residue_types: Integer that determines the number of unique residue types.
    This will determine the shape of the output profile.

Returns:
  A float Numpy array of shape (num_res, num_residue_types) with residue
  frequency (residue type count normalized by MSA depth) for every column of the
  MSA.
)";

}  // namespace

namespace alphafold3 {

void RegisterModuleMsaProfile(pybind11::module m) {
  m.def("compute_msa_profile", &ComputeMsaProfile, py::arg("msa"),
        py::arg("num_residue_types"), py::doc(kComputeMsaProfileDoc + 1));
}

}  // namespace alphafold3
