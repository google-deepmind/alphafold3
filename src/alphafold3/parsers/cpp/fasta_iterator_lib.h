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

// A C++ implementation of a FASTA parser.
#ifndef ALPHAFOLD3_SRC_ALPHAFOLD3_PARSERS_PYTHON_FASTA_ITERATOR_LIB_H_
#define ALPHAFOLD3_SRC_ALPHAFOLD3_PARSERS_PYTHON_FASTA_ITERATOR_LIB_H_

#include <fstream>
#include <ios>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"

namespace alphafold3 {

// Parse FASTA string and return list of strings with amino acid sequences.
// Returns a list of amino acid sequences only.
std::vector<std::string> ParseFasta(absl::string_view fasta_string);

// Parse FASTA string and return list of strings with amino acid sequences.
// Returns two lists: The first one with amino acid sequences, the second with
// the descriptions associated with each sequence.
std::pair<std::vector<std::string>, std::vector<std::string>>
ParseFastaIncludeDescriptions(absl::string_view fasta_string);

// Lazy FASTA parser for memory efficient FASTA parsing from a path.
class FastaFileIterator {
 public:
  // Initialise FastaFileIterator with filename of fasta. If you initialize
  // reader_ with an invalid path or empty file, it won't fail, only
  // riegeli::ReadLine within the Next method will then return false. That will
  // then trigger the "Invalid FASTA file" error.
  explicit FastaFileIterator(absl::string_view fasta_path)
      : filename_(fasta_path),
        reader_(filename_, std::ios::in),
        has_next_(true) {}

  // Returns whether there are more sequences. Returns true before first call to
  // next even if the file is empty.
  bool HasNext() const { return has_next_; }

  // Fetches the next (sequence, description) from the file.
  absl::StatusOr<std::pair<std::string, std::string>> Next();

 private:
  // Use riegeli::FileReader instead of FileLineIterator for about 2x speedup.
  std::string filename_;
  std::fstream reader_;
  std::optional<std::string> description_;
  std::string sequence_;
  bool has_next_;
};

// Lazy FASTA parser for memory efficient FASTA parsing from a string.
class FastaStringIterator {
 public:
  // Initialise FastaStringIterator with a string_view of a FASTA. If you
  // initialize it with an invalid FASTA string, it won't fail, the Next method
  // will then return false. That will then trigger the "Invalid FASTA" error.
  // WARNING: The object backing the fasta_string string_view must not be
  // deleted while this Iterator is alive.
  explicit FastaStringIterator(absl::string_view fasta_string)
      : fasta_string_(fasta_string), has_next_(true) {}

  // Returns whether there are more sequences. Returns true before first call to
  // next even if the string is empty.
  bool HasNext() const { return has_next_; }

  // Fetches the next (sequence, description) from the string.
  absl::StatusOr<std::pair<std::string, std::string>> Next();

 private:
  absl::string_view fasta_string_;
  bool has_next_;
  std::optional<std::string> description_;
  std::string sequence_;
};

}  // namespace alphafold3

#endif  // ALPHAFOLD3_SRC_ALPHAFOLD3_PARSERS_PYTHON_FASTA_ITERATOR_LIB_H_
