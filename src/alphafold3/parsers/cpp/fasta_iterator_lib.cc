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

#include "alphafold3/parsers/cpp/fasta_iterator_lib.h"

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/ascii.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "absl/strings/strip.h"

namespace alphafold3 {

// Parse FASTA string and return list of strings with amino acid sequences.
// Returns a list of amino acid sequences only.
std::vector<std::string> ParseFasta(absl::string_view fasta_string) {
  std::vector<std::string> sequences;
  std::string* sequence = nullptr;
  for (absl::string_view line_raw : absl::StrSplit(fasta_string, '\n')) {
    absl::string_view line = absl::StripAsciiWhitespace(line_raw);
    if (absl::ConsumePrefix(&line, ">")) {
      sequence = &sequences.emplace_back();
    } else if (!line.empty() && sequence != nullptr) {
      absl::StrAppend(sequence, line);
    }
  }
  return sequences;
}

// Parse FASTA string and return list of strings with amino acid sequences.
// Returns two lists: The first one with amino acid sequences, the second with
// the descriptions associated with each sequence.
std::pair<std::vector<std::string>, std::vector<std::string>>
ParseFastaIncludeDescriptions(absl::string_view fasta_string) {
  std::pair<std::vector<std::string>, std::vector<std::string>> result;
  auto& [sequences, descriptions] = result;
  std::string* sequence = nullptr;
  for (absl::string_view line_raw : absl::StrSplit(fasta_string, '\n')) {
    absl::string_view line = absl::StripAsciiWhitespace(line_raw);
    if (absl::ConsumePrefix(&line, ">")) {
      descriptions.emplace_back(line);
      sequence = &sequences.emplace_back();
    } else if (!line.empty() && sequence != nullptr) {
      absl::StrAppend(sequence, line);
    }
  }
  return result;
}

absl::StatusOr<std::pair<std::string, std::string>> FastaFileIterator::Next() {
  std::string line_str;
  while (std::getline(reader_, line_str)) {
    absl::string_view line = line_str;
    line = absl::StripAsciiWhitespace(line);
    if (absl::ConsumePrefix(&line, ">")) {
      if (!description_.has_value()) {
        description_ = line;
      } else {
        std::pair<std::string, std::string> output(sequence_, *description_);
        description_ = line;
        sequence_ = "";
        return output;
      }
    } else if (description_.has_value()) {
      absl::StrAppend(&sequence_, line);
    }
  }
  has_next_ = false;
  reader_.close();
  if (description_.has_value()) {
    return std::pair(sequence_, *description_);
  } else {
    return absl::InvalidArgumentError(
        absl::StrCat("Invalid FASTA file: ", filename_));
  }
}

absl::StatusOr<std::pair<std::string, std::string>>
FastaStringIterator::Next() {
  size_t consumed = 0;
  for (absl::string_view line_raw : absl::StrSplit(fasta_string_, '\n')) {
    consumed += line_raw.size() + 1;  // +1 for the newline character.
    absl::string_view line = absl::StripAsciiWhitespace(line_raw);
    if (absl::ConsumePrefix(&line, ">")) {
      if (!description_.has_value()) {
        description_ = line;
      } else {
        std::pair<std::string, std::string> output(sequence_, *description_);
        description_ = line;
        sequence_ = "";
        fasta_string_.remove_prefix(consumed);
        return output;
      }
    } else if (description_.has_value()) {
      absl::StrAppend(&sequence_, line);
    }
  }
  has_next_ = false;
  if (description_.has_value()) {
    return std::pair(sequence_, *description_);
  } else {
    return absl::InvalidArgumentError("Invalid FASTA string");
  }
}

}  // namespace alphafold3
