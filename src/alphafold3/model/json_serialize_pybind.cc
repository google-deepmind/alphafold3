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
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "pybind11/cast.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11_abseil/absl_casters.h"

namespace alphafold3 {
namespace {

namespace py = pybind11;

constexpr std::array<std::array<char, 5>, 10'000> MakeStringFrom2DigitNum() {
  std::array<std::array<char, 5>, 10'000> decimals = {};
  for (int i = 0; i < decimals.size(); ++i) {
    decimals[i][0] = (i / 1000 > 0) ? i / 1000 + '0' : ' ';
    decimals[i][1] = (i / 100 % 10) + '0';
    decimals[i][2] = '.';
    decimals[i][3] = (i / 10 % 10) + '0';
    decimals[i][4] = (i % 10) + '0';
  }
  return decimals;
}

static constexpr auto kStringFrom2DigitNum = MakeStringFrom2DigitNum();

// Converts floats in the range [0, 1] to strings with 2 decimal places.
// Returns "null" for NaN otherwise clamps to [0.00, 1.00].
constexpr absl::string_view StringFromFraction(float value) {
  if (std::isnan(value)) {
    return "null";
  }
  float value_clamped = std::clamp(value, 0.f, 1.f);
  const int index = static_cast<int>(std::round(value_clamped * 100.0f));
  const auto& result_array = kStringFrom2DigitNum[index];
  absl::string_view result(result_array.data() + 1, result_array.size() - 1);
  if (result.back() == '0') {
    result.remove_suffix(1);
  }
  return result;
}

// Converts floats in the range [0, 99.9] to to strings with 1 decimal places.
// Returns "null" for NaN otherwise clamps to [0, 99.9].
constexpr absl::string_view StringFrom2DigitNum1Decimal(float value) {
  if (std::isnan(value)) {
    return "null";
  }
  float value_clamped = std::clamp(value, 0.f, 99.9f);
  const int index = static_cast<int>(std::round(value_clamped * 10.0f)) * 10;
  const auto& result_array = kStringFrom2DigitNum[index];
  absl::string_view result(result_array.data(), result_array.size() - 1);
  if (result[0] == ' ') {
    result.remove_prefix(1);
  }
  return result;
}

// Converts floats in the range [0, 99.9] to to strings with 2 decimal places.
// Returns "null" for NaN otherwise clamps to [0, 99.99].
constexpr absl::string_view StringFrom2DigitNum2Decimal(float value) {
  if (std::isnan(value)) {
    return "null";
  }
  float value_clamped = std::clamp(value, 0.f, 99.99f);
  const int index = static_cast<int>(std::round(value_clamped * 100.0f));
  const auto& result_array = kStringFrom2DigitNum[index];
  absl::string_view result(result_array.data(), result_array.size());
  if (result[0] == ' ') {
    result.remove_prefix(1);
  }
  if (result.back() == '0') {
    result.remove_suffix(1);
  }
  return result;
}

constexpr absl::string_view kAtomChainIdsJson = R"(  "atom_chain_ids": [")";
constexpr absl::string_view kAtomPlddtsJson = R"(  "atom_plddts": [)";
constexpr absl::string_view kContactProbsJson = R"(  "contact_probs": [)";
constexpr absl::string_view kPaeJson = R"(  "pae": [)";
constexpr absl::string_view kTokenChainIdsJson = R"(  "token_chain_ids": [")";
constexpr absl::string_view kTokenResIdsJson = R"(  "token_res_ids": [)";

using NpFloatArray =
    py::array_t<const float, py::array::forcecast | py::array::c_style>;
using NpIntArray =
    py::array_t<const int, py::array::forcecast | py::array::c_style>;

std::string StructureConfidenceFullToJson(
    NpFloatArray pae_array,                                 //
    const std::vector<absl::string_view>& token_chain_ids,  //
    NpIntArray token_res_ids_array,                         //
    NpFloatArray atom_plddts_array,                         //
    const std::vector<absl::string_view>& atom_chain_ids,   //
    NpFloatArray contact_probs_array) {
  if (pae_array.ndim() != 2) {
    throw py::value_error("pae_array must be 2D");
  }

  if (contact_probs_array.ndim() != 2) {
    throw py::value_error("contact_probs_array must be 2D");
  }

  const size_t pae_num_rows = pae_array.shape()[0];
  const size_t pae_num_cols = pae_array.shape()[1];
  const size_t contact_probs_num_rows = contact_probs_array.shape()[0];
  const size_t contact_probs_num_cols = contact_probs_array.shape()[1];

  absl::Span<const float> pae(pae_array);
  absl::Span<const int> token_res_ids(token_res_ids_array);
  absl::Span<const float> atom_plddts(atom_plddts_array);
  absl::Span<const float> contact_probs(contact_probs_array);

  py::gil_scoped_release gil_release;

  // \b deletes the previous character comma.
  size_t capacity = 2;                       //  {\n
  capacity += kAtomChainIdsJson.size() + 1;  // [
  for (const auto& s : atom_chain_ids) {
    capacity += s.size() + 3;  // '","'
  }
  // \b\b],

  capacity += kAtomPlddtsJson.size() + 1;  // [
  capacity += atom_plddts.size() * 6;      // 10.00,
  capacity += 1;                           // \b],

  capacity += kContactProbsJson.size() + 1;  // [
  capacity += contact_probs.size() * 5;      // 0.05,
  capacity += 2 * contact_probs_num_rows;    // \b],[
  capacity += 1;                             // \b\b]],

  capacity += kPaeJson.size() + 2;  // [[
  capacity += pae.size() * 5;       // 10.0,
  capacity += 2 * pae_num_rows;     // \b],[
  // \b\b]]

  capacity += kTokenChainIdsJson.size() + 1;  // [
  for (const auto& s : token_chain_ids) {
    capacity += s.size() + 3;  // '","'
  }
  // \b]

  capacity += kTokenResIdsJson.size() + 1;  // [
  capacity += token_res_ids.size() * 2;     // 0,
  for (int res_id : token_res_ids) {
    while ((res_id /= 10) > 0) {
      capacity += 1;  // Bump for each extra digit.
    }
  }
  // \b]

  capacity += 1;  // "}"

  std::string json;
  json.reserve(capacity);
  json.append("{\n");

  // atom_chain_ids
  json.append(kAtomChainIdsJson);
  for (size_t i = 0; i + 1 < atom_chain_ids.size(); ++i) {
    json.append(atom_chain_ids[i]);
    json.append("\",\"");
  }
  json.append(atom_chain_ids.back());
  json.append("\"],\n");

  // atom_plddts
  json.append(kAtomPlddtsJson);
  for (size_t i = 0; i + 1 < atom_plddts.size(); ++i) {
    json.append(StringFrom2DigitNum2Decimal(atom_plddts[i]));
    json.append(",");
  }
  json.append(StringFrom2DigitNum2Decimal(atom_plddts.back()));
  json.append("],\n");

  // contact_probs
  json.append(kContactProbsJson);
  size_t pos = 0;
  for (size_t row = 0; row < contact_probs_num_rows; ++row) {
    json.append("[");
    for (size_t col = 0; col < contact_probs_num_cols - 1; ++col) {
      json.append(StringFromFraction(contact_probs[pos++]));
      json.append(",");
    }
    json.append(StringFromFraction(contact_probs[pos++]));
    json.append(row != contact_probs_num_rows - 1 ? "]," : "]");
  }
  json.append("],\n");

  // pae
  json.append(kPaeJson);
  pos = 0;
  for (size_t row = 0; row < pae_num_rows; ++row) {
    json.append("[");
    for (size_t col = 0; col + 1 < pae_num_cols; ++col) {
      json.append(StringFrom2DigitNum1Decimal(pae[pos++]));
      json.append(",");
    }
    json.append(StringFrom2DigitNum1Decimal(pae[pos++]));
    json.append(row != pae_num_rows - 1 ? "]," : "]");
  }
  json.append("],\n");

  // token_chain_ids
  json.append(kTokenChainIdsJson);
  for (size_t col = 0; col + 1 < token_chain_ids.size(); ++col) {
    json.append(token_chain_ids[col]);
    json.append("\",\"");
  }
  json.append(token_chain_ids.back());
  json.append("\"],\n");

  // token_res_ids
  json.append(kTokenResIdsJson);
  for (size_t col = 0; col + 1 < token_res_ids.size(); ++col) {
    json.append(absl::StrCat(token_res_ids[col]));
    json.append(",");
  }
  json.append(absl::StrCat(token_res_ids.back()));
  json.append("]\n");

  json.append("}");
  return json;
}

constexpr char kStructureConfidenceFullToJson[] = R"(
  Serializes StructureConfidenceFull to JSON faster than Python (yes, low bar).
  )";
}  // namespace

void RegisterModuleJsonSerialize(pybind11::module m) {
  m.def("structure_confidence_full_to_json",         //
        &StructureConfidenceFullToJson,              //
        py::arg("pae"),                              //
        py::arg("token_chain_ids"),                  //
        py::arg("token_res_ids"),                    //
        py::arg("atom_plddts"),                      //
        py::arg("atom_chain_ids"),                   //
        py::arg("contact_probs"),                    //
        py::doc(kStructureConfidenceFullToJson + 1)  //
  );
}

}  // namespace alphafold3
