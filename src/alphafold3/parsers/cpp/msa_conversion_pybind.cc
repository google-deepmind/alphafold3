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
#include <stdexcept>
#include <string>
#include <vector>

#include "absl/strings/ascii.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace {

namespace py = pybind11;

std::vector<std::string> ConvertA3MToStockholm(
    std::vector<absl::string_view> a3m_sequences) {
  std::vector<std::string> stockholm_sequences(a3m_sequences.size());
  auto max_length_element =
      std::max_element(a3m_sequences.begin(), a3m_sequences.end(),
                       [](absl::string_view lhs, absl::string_view rhs) {
                         return lhs.size() < rhs.size();
                       });

  for (auto& out : stockholm_sequences) {
    out.reserve(max_length_element->size());
  }

  // While any sequence has remaining columns.
  while (std::any_of(a3m_sequences.begin(), a3m_sequences.end(),
                     [](absl::string_view in) { return !in.empty(); })) {
    if (std::any_of(a3m_sequences.begin(), a3m_sequences.end(),
                    [](absl::string_view in) {
                      return !in.empty() && absl::ascii_islower(in.front());
                    })) {
      // Insertion(s) found at column.
      for (std::size_t i = 0; i < a3m_sequences.size(); ++i) {
        absl::string_view& in = a3m_sequences[i];
        std::string& out = stockholm_sequences[i];
        if (!in.empty() && absl::ascii_islower(in.front())) {
          // Consume insertion.
          out.push_back(absl::ascii_toupper(in.front()));
          in.remove_prefix(1);
        } else {
          // Row requires padding.
          out.push_back('-');
        }
      }
    } else {
      // No insertions found.
      for (std::size_t i = 0; i < a3m_sequences.size(); ++i) {
        absl::string_view& in = a3m_sequences[i];
        std::string& out = stockholm_sequences[i];
        if (!in.empty()) {
          // Consume entire column.
          out.push_back(in.front());
          in.remove_prefix(1);
        } else {
          // One alignment is shorter than the others. Should not happen with
          // valid A3M input.
          throw std::invalid_argument(absl::StrFormat(
              "a3m rows have inconsistent lengths; row %d has no columns left "
              "but not all rows are exhausted",
              i));
        }
      }
    }
  }
  return stockholm_sequences;
}

std::string AlignSequenceToGaplessQuery(absl::string_view sequence,
                                        absl::string_view query_sequence) {
  if (sequence.size() != query_sequence.size()) {
    throw py::value_error(
        absl::StrFormat("The sequence (%d) and the query sequence (%d) don't "
                        "have the same length.",
                        sequence.size(), query_sequence.size()));
  }
  std::string output;
  for (std::size_t residue_index = 0, sequence_length = sequence.size();
       residue_index < sequence_length; ++residue_index) {
    const char query_residue = query_sequence[residue_index];
    const char residue = sequence[residue_index];
    if (query_residue != '-') {
      // No gap in the query, so the residue is aligned.
      output += residue;
    } else if (residue == '-') {
      // Gap in both sequence and query, simply skip.
      continue;
    } else {
      // Gap only in the query, so this must be an inserted residue.
      output += absl::ascii_tolower(residue);
    }
  }
  return output;
}

constexpr char kConvertA3mToStockholm[] = R"(
Converts a list of sequences in a3m format to stockholm format sequences.

As an example if the input is:
abCD
CgD
fCDa

Then the output will be:
ABC-D-
--CGD-
F-C-DA

Args:
  a3m_sequences: A list of strings in a3m format.

Returns
  A list of strings converted to stockholm format.
)";

constexpr char kAlignSequenceToGaplessQuery[] = R"(
Aligns a sequence to a gapless query sequence.

This is useful when converting Stockholm MSA to A3M MSA. Example:
Seq  : AB--E
Query: A--DE
Output: Ab-E.

Args:
  sequence: A string containing to be aligned.
  query_sequence: A string containing the reference sequence to align to.

Returns
  The input sequence with gaps dropped where both the `sequence` and
  `query_sequence` have gaps, and sequence elements non-capitalized where the
  `query_sequence` has a gap, but the `sequence` does not.
)";

}  // namespace

namespace alphafold3 {

void RegisterModuleMsaConversion(pybind11::module m) {
  m.def("convert_a3m_to_stockholm", &ConvertA3MToStockholm,
        py::arg("a3m_sequences"), py::call_guard<py::gil_scoped_release>(),
        py::doc(kConvertA3mToStockholm + 1));
  m.def("align_sequence_to_gapless_query", &AlignSequenceToGaplessQuery,
        py::arg("sequence"), py::arg("query_sequence"),
        py::call_guard<py::gil_scoped_release>(),
        py::doc(kAlignSequenceToGaplessQuery + 1));
}

}  // namespace alphafold3
