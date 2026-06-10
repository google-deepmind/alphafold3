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

// A C++ implementation of a CIF parser. For the format specification see
// https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax
#ifndef ALPHAFOLD3_SRC_ALPHAFOLD3_PARSERS_PYTHON_CIF_DICT_LIB_H_
#define ALPHAFOLD3_SRC_ALPHAFOLD3_PARSERS_PYTHON_CIF_DICT_LIB_H_

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/node_hash_map.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"

namespace alphafold3 {

class CifDict {
 public:
  // Use absl::node_hash_map since it guarantees pointer stability.
  using Dict = absl::node_hash_map<std::string, std::vector<std::string>>;

  CifDict() = default;

  explicit CifDict(Dict dict)
      : dict_(std::make_shared<const Dict>(std::move(dict))) {}

  // Converts a CIF string into a dictionary mapping each CIF field to a list of
  // values that field contains.
  static absl::StatusOr<CifDict> FromString(absl::string_view cif_string);

  // Converts the CIF into into a string that is a valid CIF file.
  absl::StatusOr<std::string> ToString() const;

  // Extracts loop associated with a prefix from mmCIF data as a list.
  // Reference for loop_ in mmCIF:
  // http://mmcif.wwpdb.org/docs/tutorials/mechanics/pdbx-mmcif-syntax.html
  // Args:
  // prefix: Prefix shared by each of the data items in the loop.
  //   e.g. '_entity_poly_seq.', where the data items are _entity_poly_seq.num,
  //   _entity_poly_seq.mon_id. Should include the trailing period.
  //
  // Returns a list of dicts; each dict represents 1 entry from an mmCIF loop.
  // Lifetime of string_views tied to this.
  absl::StatusOr<
      std::vector<absl::flat_hash_map<absl::string_view, absl::string_view>>>
  ExtractLoopAsList(absl::string_view prefix) const;

  // Extracts loop associated with a prefix from mmCIF data as a dictionary.
  // Args:
  // prefix: Prefix shared by each of the data items in the loop.
  //   e.g. '_entity_poly_seq.', where the data items are _entity_poly_seq.num,
  //   _entity_poly_seq.mon_id. Should include the trailing period.
  // index: Which item of loop data should serve as the key.
  //
  // Returns a dict of dicts; each dict represents 1 entry from an mmCIF loop,
  // indexed by the index column.
  // Lifetime of string_views tied to this.
  absl::StatusOr<absl::flat_hash_map<
      absl::string_view,
      absl::flat_hash_map<absl::string_view, absl::string_view>>>
  ExtractLoopAsDict(absl::string_view prefix, absl::string_view index) const;

  // Returns value at key if present or an empty list.
  absl::Span<const std::string> operator[](absl::string_view key) const {
    auto it = dict_->find(key);
    if (it != dict_->end()) {
      return it->second;
    }
    return {};
  }

  // Returns boolean of whether dict contains key.
  bool Contains(absl::string_view key) const { return dict_->contains(key); }

  // Returns number of values for the given key if present, 0 otherwise.
  size_t ValueLength(absl::string_view key) const {
    return (*this)[key].size();
  }

  // Returns the size of the underlying dictionary.
  std::size_t Length() { return dict_->size(); }

  // Creates a copy of this CifDict object that will contain the original values
  // but only if not updated by the given dictionary.
  // E.g. if the CifDict = {a: [a1, a2], b: [b1]} and other = {a: [x], c: [z]},
  // you will get {a: [x], b: [b1], c: [z]}.
  CifDict CopyAndUpdate(Dict other) const {
    other.insert(dict_->begin(), dict_->end());
    return CifDict(std::move(other));
  }

  // Returns the value of the special CIF data_ field.
  absl::string_view GetDataName() const {
    // The data_ element has to be present by construction.
    if (auto it = dict_->find("data_");
        it != dict_->end() && !it->second.empty()) {
      return it->second.front();
    } else {
      return "";
    }
  }

  const std::shared_ptr<const Dict>& dict() const { return dict_; }

 private:
  std::shared_ptr<const Dict> dict_;
};

// Tokenizes a CIF string into a list of string tokens. This is more involved
// than just a simple split on whitespace as CIF allows comments and quoting.
absl::StatusOr<std::vector<std::string>> Tokenize(absl::string_view cif_string);

// Tokenizes a single line of a CIF string.
absl::StatusOr<std::vector<absl::string_view>> SplitLine(
    absl::string_view line);

// Parses a CIF string with multiple data records and returns a mapping from
// record names to CifDict objects. For instance, the following CIF string:
//
// data_001
// _foo bar
//
// data_002
// _foo baz
//
// will be parsed as:
// {'001': CifDict({'_foo': ['bar']}),
//  '002': CifDict({'_foo': ['baz']})}
absl::StatusOr<absl::flat_hash_map<std::string, CifDict>> ParseMultiDataCifDict(
    absl::string_view cif_string);

}  // namespace alphafold3

#endif  // ALPHAFOLD3_SRC_ALPHAFOLD3_PARSERS_PYTHON_CIF_DICT_LIB_H_
