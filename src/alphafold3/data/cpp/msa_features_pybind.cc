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

#include "alphafold3/data/cpp/msa_features_pybind.h"

#include <algorithm>
#include <array>
#include <string>
#include <utility>
#include <vector>

#include "absl/strings/ascii.h"
#include "absl/strings/match.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "pybind11/cast.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11_abseil/absl_casters.h"

namespace {

namespace py = pybind11;

typedef std::pair<py::array_t<int>, py::array_t<int>> MsaInfo;

constexpr std::array<int, 256> ProteinToId() {
  std::array<int, 256> char_map = {};
  for (auto& c : char_map) {
    c = -1;
  }

  char_map[static_cast<unsigned char>('A')] = 0;
  char_map[static_cast<unsigned char>('B')] = 3;  // Same as D.
  char_map[static_cast<unsigned char>('C')] = 4;
  char_map[static_cast<unsigned char>('D')] = 3;
  char_map[static_cast<unsigned char>('E')] = 6;
  char_map[static_cast<unsigned char>('F')] = 13;
  char_map[static_cast<unsigned char>('G')] = 7;
  char_map[static_cast<unsigned char>('H')] = 8;
  char_map[static_cast<unsigned char>('I')] = 9;
  char_map[static_cast<unsigned char>('J')] = 20;  // Same as unknown (X).
  char_map[static_cast<unsigned char>('K')] = 11;
  char_map[static_cast<unsigned char>('L')] = 10;
  char_map[static_cast<unsigned char>('M')] = 12;
  char_map[static_cast<unsigned char>('N')] = 2;
  char_map[static_cast<unsigned char>('O')] = 20;  // Same as unknown (X).
  char_map[static_cast<unsigned char>('P')] = 14;
  char_map[static_cast<unsigned char>('Q')] = 5;
  char_map[static_cast<unsigned char>('R')] = 1;
  char_map[static_cast<unsigned char>('S')] = 15;
  char_map[static_cast<unsigned char>('T')] = 16;
  char_map[static_cast<unsigned char>('U')] = 4;  // Same as C.
  char_map[static_cast<unsigned char>('V')] = 19;
  char_map[static_cast<unsigned char>('W')] = 17;
  char_map[static_cast<unsigned char>('X')] = 20;
  char_map[static_cast<unsigned char>('Y')] = 18;
  char_map[static_cast<unsigned char>('Z')] = 6;  // Same as E.
  char_map[static_cast<unsigned char>('-')] = 21;
  return char_map;
}

constexpr std::array<int, 256> RnaToId() {
  std::array<int, 256> char_map = {};
  for (auto& c : char_map) {
    c = -1;
  }
  for (unsigned char c = 'A'; c != 'Z' + 1; ++c) {
    char_map[c] = 30;  // Map non-standard residues to UNK_NUCLEIC (N) -> 30.
  }
  // Continue the RNA indices from where the residue indices in ProteinToId end.
  char_map[static_cast<unsigned char>('-')] = 21;
  char_map[static_cast<unsigned char>('A')] = 22;
  char_map[static_cast<unsigned char>('G')] = 23;
  char_map[static_cast<unsigned char>('C')] = 24;
  char_map[static_cast<unsigned char>('U')] = 25;
  return char_map;
}

constexpr std::array<int, 256> DnaToId() {
  std::array<int, 256> char_map = {};
  for (auto& c : char_map) {
    c = -1;
  }
  for (unsigned char c = 'A'; c != 'Z' + 1; ++c) {
    char_map[c] = 30;  // Map non-standard residues to UNK_NUCLEIC (N) -> 30.
  }
  // Continue the DNA indices from where the residue indices in RnaToId end.
  char_map[static_cast<unsigned char>('-')] = 21;
  char_map[static_cast<unsigned char>('A')] = 26;
  char_map[static_cast<unsigned char>('G')] = 27;
  char_map[static_cast<unsigned char>('C')] = 28;
  char_map[static_cast<unsigned char>('T')] = 29;
  return char_map;
}

MsaInfo ExtractMsaFeatures(std::vector<absl::string_view> msa_sequences,
                           absl::string_view chain_poly_type) {
  static constexpr auto kProteinToId = ProteinToId();
  static constexpr auto kRnaToId = RnaToId();
  static constexpr auto kDnaToId = DnaToId();

  absl::Span<const int> char_map;
  if (chain_poly_type == "polyribonucleotide") {
    char_map = kRnaToId;
  } else if (chain_poly_type == "polydeoxyribonucleotide") {
    char_map = kDnaToId;
  } else if (chain_poly_type == "polypeptide(L)") {
    char_map = kProteinToId;
  } else {
    throw py::value_error(
        absl::StrFormat("Chain type %s invalid.", chain_poly_type));
  }
  if (msa_sequences.empty()) {
    int shape[] = {0, 0};
    py::array_t<int> empty_array(shape);
    return {empty_array, empty_array};
  }
  int num_rows = msa_sequences.size();
  int num_cols = std::count_if(
      msa_sequences[0].begin(), msa_sequences[0].end(), [&char_map](char c) {
        return char_map[static_cast<unsigned char>(c)] != -1;
      });
  py::array_t<int> msa_arr({num_rows, num_cols});
  py::array_t<int> deletions_arr({num_rows, num_cols});
  {
    py::gil_scoped_release release_gil;
    absl::Span<int> msa(msa_arr.mutable_data(), msa_arr.size());
    absl::Span<int> deletions(deletions_arr.mutable_data(),
                              deletions_arr.size());

    int problem_row = 0;
    int flatten_index = 0;
    for (const auto& msa_sequence : msa_sequences) {
      int deletion_count = 0;
      int upper_count = 0;
      int problem_col = 0;
      std::vector<std::string> problems;
      for (char current : msa_sequence) {
        int msa_id = char_map[static_cast<unsigned char>(current)];
        if (msa_id == -1) {
          if (!absl::ascii_islower(current)) {
            problems.push_back(absl::StrFormat("(%d, %d):%c", problem_row,
                                               problem_col, current));
          }
          ++deletion_count;
        } else {
          if (flatten_index < deletions.size()) {
            deletions[flatten_index] = deletion_count;
            msa[flatten_index] = msa_id;
          }
          deletion_count = 0;
          ++flatten_index;
          ++upper_count;
        }
        ++problem_col;
      }
      if (!problems.empty()) {
        throw py::value_error(
            absl::StrFormat("Unknown residues in MSA: %s. target_sequence: %s",
                            absl::StrJoin(problems, ", "), msa_sequences[0]));
      }
      if (upper_count != num_cols) {
        throw py::value_error(absl::StrFormat(
            "Invalid shape all strings must have the same number "
            "of non-lowercase characters; First string has %d "
            "non-lowercase characters but '%s' has %d. target_sequence: %s",
            num_cols, msa_sequence, upper_count, msa_sequences[0]));
      }
      ++problem_row;
    }
  }

  return {std::move(msa_arr), std::move(deletions_arr)};
}

std::vector<py::bytes> ExtractSpeciesIds(
    std::vector<absl::string_view> msa_descriptions) {
  std::vector<py::bytes> species_ids;
  species_ids.reserve(msa_descriptions.size());

  for (absl::string_view msa_description : msa_descriptions) {
    // UniProtKB SwissProt/TrEMBL dbs have the following description format:
    // `db|UniqueIdentifier|EntryName`, e.g. `sp|P0C2L1|A3X1_LOXLA` or
    // `tr|A0A146SKV9|A0A146SKV9_FUNHE`.
    // See https://www.uniprot.org/help/fasta-headers for full documentation.
    msa_description = absl::StripAsciiWhitespace(msa_description);
    absl::string_view part =
        msa_description.substr(0, msa_description.find(' '));
    if (part.empty() ||
        !(part.starts_with("sp|") || part.starts_with("tr|"))) {
      species_ids.emplace_back("");
      continue;
    }
    int count = 0;
    absl::string_view selected_identifier;
    for (absl::string_view identifier : absl::StrSplit(part, '|')) {
      ++count;
      if (count == 2 && absl::StrContains(identifier, '-')) {
        // Exclude alternative isoforms, e.g. "sp|O88974-2|SETB1_MOUSE".
        break;
      } else if (count == 3) {
        selected_identifier = identifier;
      } else if (count > 3) {
        // Too many.
        selected_identifier = "";
        break;
      }
    }
    std::pair<absl::string_view, absl::string_view> id_and_species =
        absl::StrSplit(selected_identifier, '_');
    // There is an additional range after a "/" added by Jackhmmer.
    // E.g. "sp|Q32KD2|SETB1_DROME/444-700".
    species_ids.emplace_back(
        id_and_species.second.substr(0, id_and_species.second.find('/')));
  }
  return species_ids;
}

constexpr char kExtractMsaFeatures[] = R"(
Extracts MSA features which are too expensive to compute in Python.

Example:
The input raw MSA is: `[["AAAAAA"], ["Ai-CiDiiiEFa"]]`
The output MSA will be: `[["AAAAAA"], ["A-CDEF"]]`
The deletions will be: `[[0, 0, 0, 0, 0, 0], [0, 1, 0, 1, 3, 0]]`

Args:
  msa_sequences: A list of strings, each string with one MSA sequence.
    Each string must have the same, constant number of non-lowercase
    (matching) residues.
  chain_poly_type: Either 'polypeptide(L)' (protein),
    'polyribonucleotide' (RNA), or 'polydeoxyribonucleotide' (DNA). Use
    the appropriate string constant from mmcif_names.py.

Returns:
  A tuple with:
  * MSA array of shape (num_seq, num_res) that contains only the uppercase
    characters or gaps (-) from the original MSA.
  * Deletions array of shape (num_seq, num_res) that contains the number
    of deletions (lowercase letters in the MSA) to the left from each
    non-deleted residue (uppercase letters in the MSA).

Raises:
InvalidArgumentError if any of the preconditions are not met.
)";

constexpr char kExtractSpeciesIds[] = R"(
Extracts species ID from MSA UniProtKB sequence identifiers.

Args:
  msa_descriptions: The descriptions (the FASTA/A3M comment line) for each
    of the sequences.

Returns:
  Extracted UniProtKB species IDs if there is a regex match for each
  description line, blank if the regex doesn't match. Returned as bytes to
  save conversion when turning into a feature array.
)";

}  // namespace

namespace alphafold3 {

void RegisterModuleMsaFeatures(pybind11::module m) {
  m.def("extract_msa_features", &ExtractMsaFeatures, py::arg("msa_sequences"),
        py::arg("chain_poly_type"), py::doc(kExtractMsaFeatures + 1));
  m.def("extract_species_ids", &ExtractSpeciesIds, py::arg("msa_descriptions"),
        py::doc(kExtractSpeciesIds + 1));
}

}  // namespace alphafold3
