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

#include "alphafold3/parsers/cpp/cif_dict_lib.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/algorithm/container.h"
#include "absl/base/attributes.h"
#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/node_hash_map.h"
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/ascii.h"
#include "absl/strings/match.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "absl/strings/strip.h"

namespace alphafold3 {
namespace {

bool IsQuote(const char symbol) { return symbol == '\'' || symbol == '"'; }
bool IsWhitespace(const char symbol) { return symbol == ' ' || symbol == '\t'; }

// Splits line into tokens and invokes on_token for each token.
template <typename TokenFn>
ABSL_ATTRIBUTE_ALWAYS_INLINE absl::Status SplitLineInline(
    absl::string_view line, TokenFn&& on_token) {
  // See https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax
  for (int i = 0, line_length = line.length(); i < line_length;) {
    // Skip whitespace (spaces or tabs).
    while (IsWhitespace(line[i])) {
      if (++i == line_length) {
        break;
      }
    }
    if (i == line_length) {
      break;
    }

    // Skip comments (from # until the end of the line). If # is a non-comment
    // character, it must be inside a quoted token.
    if (line[i] == '#') {
      break;
    }

    int start_index;
    int end_index;
    if (IsQuote(line[i])) {
      // Token in single or double quotes. CIF v1.1 specification considers a
      // quote to be an opening quote only if it is at the beginning of a token.
      // So e.g. A' B has tokens A' and B. Also, ""A" is a token "A.
      const char quote_char = line[i++];
      start_index = i;

      // Find matching quote. The double loop is not strictly necessary, but
      // optimises a bit better.
      while (true) {
        while (i < line_length && line[i] != quote_char) {
          ++i;
        }
        if (i == line_length) {
          // Reached the end of the line while still being inside a token.
          return absl::InvalidArgumentError(
              absl::StrCat("Line ended with quote open: ", line));
        }
        if (i + 1 == line_length || IsWhitespace(line[i + 1])) {
          break;
        }
        ++i;
      }
      end_index = i++;
    } else {
      // Non-quoted token. Read until reaching whitespace.
      start_index = i++;
      while (i < line_length && !IsWhitespace(line[i])) {
        ++i;
      }
      end_index = i;
    }

    absl::Status on_token_status =
        on_token(line.substr(start_index, end_index - start_index));
    if (!on_token_status.ok()) {
      return on_token_status;
    }
  }

  return absl::OkStatus();
}

template <typename TokenFn>
ABSL_ATTRIBUTE_ALWAYS_INLINE absl::Status TokenizeInternal(
    absl::string_view cif_string, TokenFn&& on_token) {
  const std::vector<absl::string_view> lines = absl::StrSplit(cif_string, '\n');
  int line_num = 0;
  while (line_num < lines.size()) {
    auto line = absl::StripSuffix(lines[line_num], "\r");
    line_num++;

    if (line.empty() || line[0] == '#') {
      // Skip empty lines or lines that contain only comments.
      continue;
    } else if (line[0] == ';') {
      // Leading whitespace on each line must be preserved while trailing
      // whitespace may be stripped.
      std::vector<absl::string_view> multiline_tokens;
      // Strip the leading ";".
      multiline_tokens.push_back(
          absl::StripTrailingAsciiWhitespace(line.substr(1)));
      while (line_num < lines.size()) {
        auto multiline = absl::StripTrailingAsciiWhitespace(lines[line_num]);
        line_num++;
        if (!multiline.empty() && multiline[0] == ';') {
          break;
        } else if (line_num == lines.size()) {
          return absl::InvalidArgumentError(
              "Last multiline token is not terminated by a semicolon.");
        }
        multiline_tokens.push_back(multiline);
      }
      absl::Status on_token_status =
          on_token(absl::StrJoin(multiline_tokens, "\n"));
      if (!on_token_status.ok()) {
        return on_token_status;
      }
    } else {
      absl::Status on_token_status = SplitLineInline(line, on_token);
      if (!on_token_status.ok()) {
        return on_token_status;
      }
    }
  }
  return absl::OkStatus();
}

// Returns whether the token doesn't need any quoting. This is true if the token
// isn't empty and contains only safe characters [A-Za-z0-9.?-].
bool IsTrivialToken(const absl::string_view value) {
  if (value.empty()) {
    return false;
  }

  return absl::c_all_of(value, [](char c) {
    return absl::ascii_isalnum(c) || c == '.' || c == '?' || c == '-';
  });
}

// Returns whether the token needs to be a multiline token. This happens if it
// has a newline or both single and double quotes.
bool IsMultiLineToken(const absl::string_view value) {
  bool has_single_quotes = false;
  bool has_double_quotes = false;
  for (const char c : value) {
    if (c == '\n') {
      return true;
    } else if (c == '\'') {
      has_single_quotes = true;
    } else if (c == '"') {
      has_double_quotes = true;
    }
  }
  return has_single_quotes && has_double_quotes;
}

absl::string_view GetEscapeQuote(const absl::string_view value) {
  // Empty values should not happen, but if so, they should be quoted.
  if (value.empty()) {
    return "\"";
  }

  // The value must not start with one of these CIF keywords.
  if (absl::StartsWithIgnoreCase(value, "data_") ||
      absl::StartsWithIgnoreCase(value, "loop_") ||
      absl::StartsWithIgnoreCase(value, "save_") ||
      absl::StartsWithIgnoreCase(value, "stop_") ||
      absl::StartsWithIgnoreCase(value, "global_")) {
    return "\"";
  }

  // The first character must not be a special character.
  const char first = value.front();
  if (first == '_' || first == '#' || first == '$' || first == '[' ||
      first == ']' || first == ';') {
    return "\"";
  }

  // No quotes or whitespace allowed inside. Rare case when both double and
  // single quotes are present is handled by IsMultiLineToken.
  bool use_double_quote = true;
  bool use_single_quote = true;
  bool needs_quote = false;
  for (const char c : value) {
    if (c == ' ' || c == '\t') {
      needs_quote = true;
    } else if (c == '"') {
      needs_quote = true;
      use_double_quote = false;
    } else if (c == '\'') {
      needs_quote = true;
      use_single_quote = false;
    }
  }
  if (needs_quote && use_double_quote) {
    return "\"";
  } else if (needs_quote && use_single_quote) {
    return "'";
  }
  return "";
}

int RecordIndex(absl::string_view record) {
  if (record == "_entry") {
    return 0;  // _entry is always first.
  }
  if (record == "_atom_site") {
    return 2;  // _atom_site is always last.
  }
  return 1;  // other records are between _entry and _atom_site.
}

struct RecordOrder {
  using is_transparent = void;  // Enable heterogeneous lookup.
  bool operator()(absl::string_view lhs, absl::string_view rhs) const {
    std::size_t lhs_index = RecordIndex(lhs);
    std::size_t rhs_index = RecordIndex(rhs);
    return std::tie(lhs_index, lhs) < std::tie(rhs_index, rhs);
  }
};

// Make sure the _atom_site loop columns are sorted in the PDB-standard way.
constexpr absl::string_view kAtomSiteSortOrder[] = {
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_entity_id",
    "_atom_site.label_seq_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.pdbx_formal_charge",
    "_atom_site.auth_seq_id",
    "_atom_site.auth_comp_id",
    "_atom_site.auth_asym_id",
    "_atom_site.auth_atom_id",
    "_atom_site.pdbx_PDB_model_num",
};

size_t AtomSiteIndex(absl::string_view atom_site) {
  return std::distance(std::begin(kAtomSiteSortOrder),
                       absl::c_find(kAtomSiteSortOrder, atom_site));
}

struct AtomSiteOrder {
  bool operator()(absl::string_view lhs, absl::string_view rhs) const {
    auto lhs_index = AtomSiteIndex(lhs);
    auto rhs_index = AtomSiteIndex(rhs);
    return std::tie(lhs_index, lhs) < std::tie(rhs_index, rhs);
  }
};

class Column {
 public:
  Column(absl::string_view key, const std::vector<std::string>* values)
      : key_(key), values_(values) {
    int max_value_length = 0;
    for (size_t i = 0; i < values->size(); ++i) {
      absl::string_view value = (*values)[i];
      if (IsTrivialToken(value)) {
        // Shortcut for the most common cases where no quoting/multiline needed.
        max_value_length = std::max<int>(max_value_length, value.size());
        continue;
      } else if (IsMultiLineToken(value)) {
        values_with_newlines_.insert(i);
      } else {
        absl::string_view quote = GetEscapeQuote(value);
        if (!quote.empty()) {
          values_with_quotes_[i] = quote;
        }
        max_value_length =
            std::max<int>(max_value_length, value.size() + quote.size() * 2);
      }
    }
    max_value_length_ = max_value_length;
  }

  absl::string_view key() const { return key_; }

  const std::vector<std::string>* values() const { return values_; }

  int max_value_length() const { return max_value_length_; }

  bool has_newlines(size_t index) const {
    return values_with_newlines_.contains(index);
  }

  absl::string_view quote(size_t index) const {
    if (auto it = values_with_quotes_.find(index);
        it != values_with_quotes_.end()) {
      return it->second;
    }
    return "";
  }

 private:
  absl::string_view key_;
  const std::vector<std::string>* values_;
  int max_value_length_;
  // Values with newlines or quotes are very rare in a typical CIF file.
  absl::flat_hash_set<size_t> values_with_newlines_;
  absl::flat_hash_map<size_t, absl::string_view> values_with_quotes_;
};

struct GroupedKeys {
  std::vector<Column> grouped_columns;
  int max_key_length;
  int value_size;
};

absl::Status CheckLoopColumnSizes(int num_loop_keys, int num_loop_values) {
  if ((num_loop_keys > 0) && (num_loop_values % num_loop_keys != 0)) {
    return absl::InvalidArgumentError(absl::StrFormat(
        "The number of values (%d) in a loop is not a multiple of the "
        "number of the loop's columns (%d)",
        num_loop_values, num_loop_keys));
  }
  return absl::OkStatus();
}

}  // namespace

absl::StatusOr<CifDict> CifDict::FromString(absl::string_view cif_string) {
  CifDict::Dict cif;

  bool is_first_token = true;
  bool loop_flag = false;
  std::string key;

  // Counters for CIF loop_ regions.
  int num_loop_keys = 0;
  // Loops have usually O(10) columns but could have up to O(10^6) rows. It is
  // therefore wasteful to look up the cif vector where to add a loop value
  // since that means doing `columns * rows` map lookups. If we save pointers to
  // these loop column fields instead, we need only 1 cif lookup per column.
  std::vector<std::vector<std::string>*> loop_column_values;
  // Use a resetting column index instead of computing modulo on every token.
  // This is the hot path for _atom_site tables with millions of values.
  int token_column_index = 0;
  // Total number of loop values seen (for column size validation).
  int loop_token_count = 0;

  auto process_token = [&](absl::string_view token)
                           ABSL_ATTRIBUTE_ALWAYS_INLINE -> absl::Status {
    if (is_first_token) {
      is_first_token = false;
      // The first token should be data_XXX. Split into key = data, value = XXX.
      if (!absl::ConsumePrefix(&token, "data_")) {
        return absl::InvalidArgumentError(
            "The CIF file does not start with the data_ field.");
      }
      if (token.empty()) {
        return absl::InvalidArgumentError(
            "The CIF file does not contain a data block name.");
      }
      cif["data_"].emplace_back(token);
      return absl::OkStatus();
    }
    if ((token.size() == 5 && (token[0] == 'l' || token[0] == 'L')) &&
        absl::EqualsIgnoreCase(token, "loop_")) {
      // A new loop started, check the previous loop and get rid of its data.
      absl::Status loop_status =
          CheckLoopColumnSizes(num_loop_keys, loop_token_count);
      if (!loop_status.ok()) {
        return loop_status;
      }
      loop_flag = true;
      loop_column_values.clear();
      token_column_index = 0;
      loop_token_count = 0;
      num_loop_keys = 0;
      return absl::OkStatus();
    } else if (loop_flag) {
      // The second condition checks we are in the first column. Some mmCIF
      // files (e.g. 4q9r) have values in later columns starting with an
      // underscore and we don't want to read these as keys.
      if (token_column_index == 0 && !token.empty() && token[0] == '_') {
        if (loop_token_count > 0) {
          // We are out of the loop.
          loop_flag = false;
        } else {
          // We are in the keys (column names) section of the loop.
          auto [it, inserted] = cif.try_emplace(token);
          if (!inserted) {
            return absl::InvalidArgumentError(
                absl::StrCat("Duplicate loop key: '", token, "'"));
          }
          auto& columns = it->second;
          columns.clear();

          // Heuristic: _atom_site is typically the largest table in an mmCIF
          // with ~16 columns. Make sure we reserve enough space for its values.
          if (absl::StartsWith(token, "_atom_site.")) {
            columns.reserve(cif_string.size() / 80);
          }

          // Save the pointer to the loop column values.
          loop_column_values.push_back(&columns);
          num_loop_keys += 1;
          return absl::OkStatus();
        }
      } else {
        // We are in the values section of the loop. We have a pointer to the
        // loops' values, add the new token in there.
        if (token_column_index >= loop_column_values.size()) {
          return absl::InvalidArgumentError(
              absl::StrCat("Too many columns at: '", token,
                           "' at column index: ", token_column_index,
                           " expected at most: ", loop_column_values.size()));
        }
        loop_column_values[token_column_index]->emplace_back(token);
        loop_token_count++;
        if (++token_column_index == num_loop_keys) {
          // We have completed a row, reset the column index for the next row.
          token_column_index = 0;
        }
        return absl::OkStatus();
      }
    }
    if (key.empty()) {
      key = std::string(token);
      if (!absl::StartsWith(key, "_")) {
        return absl::InvalidArgumentError(
            absl::StrCat("Key '", key, "' does not start with an underscore."));
      }
    } else {
      auto [it, inserted] = cif.try_emplace(key);
      if (!inserted) {
        return absl::InvalidArgumentError(
            absl::StrCat("Duplicate key: '", key, "'"));
      }
      it->second.emplace_back(token);
      key = "";
    }
    return absl::OkStatus();
  };
  absl::Status tokenize_status = TokenizeInternal(cif_string, process_token);
  if (!tokenize_status.ok()) {
    return tokenize_status;
  }
  if (is_first_token) {
    return absl::InvalidArgumentError("The CIF file must not be empty.");
  }
  absl::Status loop_status =
      CheckLoopColumnSizes(num_loop_keys, loop_token_count);
  if (!loop_status.ok()) {
    return loop_status;
  }
  return CifDict(std::move(cif));
}

absl::StatusOr<std::string> CifDict::ToString() const {
  std::string output;

  absl::string_view data_name;
  // Check that the data_ field exists.
  if (auto name_it = (*dict_).find("data_");
      name_it == (*dict_).end() || name_it->second.empty()) {
    return absl::InvalidArgumentError(
        "The CIF must contain a valid name for this data block in the special "
        "data_ field.");
  } else {
    data_name = name_it->second.front();
  }

  if (absl::c_any_of(data_name,
                     [](char i) { return absl::ascii_isspace(i); })) {
    return absl::InvalidArgumentError(absl::StrFormat(
        "The CIF data block name must not contain any whitespace characters, "
        "got '%s'.",
        data_name));
  }
  absl::StrAppend(&output, "data_", data_name, "\n#\n");

  // Group keys by their prefix. Use btree_map to iterate in alphabetical order,
  // but with some keys being placed at the end (e.g. _atom_site).
  absl::btree_map<std::string, GroupedKeys, RecordOrder> grouped_keys;
  for (const auto& [key, values] : *dict_) {
    if (key == "data_") {
      continue;  // Skip the special data_ key, we are already done with it.
    }
    const std::pair<absl::string_view, absl::string_view> key_parts =
        absl::StrSplit(key, absl::MaxSplits('.', 1));
    const absl::string_view key_prefix = key_parts.first;
    auto [it, inserted] = grouped_keys.emplace(key_prefix, GroupedKeys{});
    GroupedKeys& grouped_key = it->second;
    grouped_key.grouped_columns.push_back(Column(key, &values));
    if (inserted) {
      grouped_key.max_key_length = key.length();
      grouped_key.value_size = values.size();
    } else {
      grouped_key.max_key_length =
          std::max<int>(key.length(), grouped_key.max_key_length);
      if (grouped_key.value_size != values.size()) {
        return absl::InvalidArgumentError(
            absl::StrFormat("Values for key %s have different length (%d) than "
                            "the other values with the same key prefix (%d).",
                            key, values.size(), grouped_key.value_size));
      }
    }
  }

  for (auto& [key_prefix, group_info] : grouped_keys) {
    if (key_prefix == "_atom_site") {
      // Make sure we sort the _atom_site loop in the standard way.
      absl::c_sort(group_info.grouped_columns,
                   [](const Column& lhs, const Column& rhs) {
                     return AtomSiteOrder{}(lhs.key(), rhs.key());
                   });
    } else {
      // Make the key ordering within a key group deterministic.
      absl::c_sort(group_info.grouped_columns,
                   [](const Column& lhs, const Column& rhs) {
                     return lhs.key() < rhs.key();
                   });
    }

    // Force `_atom_site` field to always be a loop. This resolves issues with
    // third party mmCIF parsers such as OpenBabel which always expect a loop
    // even when there is only a single atom present.
    if (group_info.value_size == 1 && key_prefix != "_atom_site") {
      // Plain key-value pairs, output them as they are.
      for (const Column& grouped_column : group_info.grouped_columns) {
        int width = group_info.max_key_length + 1;
        size_t start_pos = output.size();
        output.append(width, ' ');
        auto out_it = output.begin() + start_pos;
        absl::c_copy(grouped_column.key(), out_it);
        // Append the value, handle multi-line/quoting.
        absl::string_view value = grouped_column.values()->front();
        if (grouped_column.has_newlines(0)) {
          absl::StrAppend(&output, "\n;", value, "\n;\n");  // Multi-line value.
        } else {
          const absl::string_view quote_char = grouped_column.quote(0);
          absl::StrAppend(&output, quote_char, value, quote_char, "\n");
        }
      }
    } else {
      // CIF loop. Output the column names, then the rows with data.
      absl::StrAppend(&output, "loop_\n");
      for (Column& grouped_column : group_info.grouped_columns) {
        absl::StrAppend(&output, grouped_column.key(), "\n");
      }
      // Write the loop values, line by line. This is the most expensive part
      // since this path is taken to write the entire atom site table which has
      // about 20 columns, but thousands of rows.
      for (int i = 0; i < group_info.value_size; i++) {
        for (int column_index = 0;
             column_index < group_info.grouped_columns.size(); ++column_index) {
          const Column& grouped_column =
              group_info.grouped_columns[column_index];
          const absl::string_view value = (*grouped_column.values())[i];
          if (grouped_column.has_newlines(i)) {
            // Multi-line. This is very rarely taken path.
            if (column_index == 0) {
              // No extra newline before leading ;, already inserted.
              absl::StrAppend(&output, ";", value, "\n;\n");
            } else if (column_index == group_info.grouped_columns.size() - 1) {
              // No extra newline after trailing ;, will be inserted.
              absl::StrAppend(&output, "\n;", value, "\n;");
            } else {
              absl::StrAppend(&output, "\n;", value, "\n;\n");
            }
          } else {
            size_t start_pos = output.size();
            output.append(grouped_column.max_value_length() + 1, ' ');
            auto out_it = output.begin() + start_pos;
            absl::string_view quote = grouped_column.quote(i);
            if (!quote.empty()) {
              out_it = absl::c_copy(quote, out_it);
              out_it = absl::c_copy(value, out_it);
              absl::c_copy(quote, out_it);
            } else {
              absl::c_copy(value, out_it);
            }
          }
        }
        absl::StrAppend(&output, "\n");
      }
    }
    absl::StrAppend(&output, "#\n");  // Comment token after every key group.
  }
  return output;
}

absl::StatusOr<
    std::vector<absl::flat_hash_map<absl::string_view, absl::string_view>>>
CifDict::ExtractLoopAsList(absl::string_view prefix) const {
  std::vector<absl::string_view> column_names;
  std::vector<std::vector<absl::string_view>> column_data;

  for (const auto& element : *dict_) {
    if (absl::StartsWith(element.first, prefix)) {
      column_names.emplace_back(element.first);
      auto& cells = column_data.emplace_back();
      cells.insert(cells.begin(), element.second.begin(), element.second.end());
    }
  }
  // Make sure all columns have the same number of rows.
  const std::size_t num_rows = column_data.empty() ? 0 : column_data[0].size();
  for (const auto& column : column_data) {
    if (column.size() != num_rows) {
      return absl::InvalidArgumentError(absl::StrCat(
          GetDataName(),
          ": Columns do not have the same number of rows for prefix: '", prefix,
          "'. One possible reason could be not including the trailing dot, "
          "e.g. '_atom_site.'."));
    }
  }

  std::vector<absl::flat_hash_map<absl::string_view, absl::string_view>> result;
  result.reserve(num_rows);
  CHECK_EQ(column_names.size(), column_data.size());
  for (std::size_t row_index = 0; row_index < num_rows; ++row_index) {
    auto& row_dict = result.emplace_back();
    row_dict.reserve(column_names.size());
    for (int col_index = 0; col_index < column_names.size(); ++col_index) {
      row_dict[column_names[col_index]] = column_data[col_index][row_index];
    }
  }
  return result;
}

absl::StatusOr<absl::flat_hash_map<
    absl::string_view,
    absl::flat_hash_map<absl::string_view, absl::string_view>>>
CifDict::ExtractLoopAsDict(absl::string_view prefix,
                           absl::string_view index) const {
  if (!absl::StartsWith(index, prefix)) {
    return absl::InvalidArgumentError(
        absl::StrCat(GetDataName(), ": The loop index '", index,
                     "' must start with the loop prefix '", prefix, "'."));
  }
  absl::flat_hash_map<absl::string_view,
                      absl::flat_hash_map<absl::string_view, absl::string_view>>
      result;
  auto loop_as_list = ExtractLoopAsList(prefix);
  if (!loop_as_list.ok()) {
    return loop_as_list.status();
  }
  result.reserve(loop_as_list->size());
  for (auto& entry : *loop_as_list) {
    if (const auto it = entry.find(index); it != entry.end()) {
      result[it->second] = entry;
    } else {
      return absl::InvalidArgumentError(absl::StrCat(
          GetDataName(), ": The index column '", index,
          "' could not be found in the loop with prefix '", prefix, "'."));
    }
  }
  return result;
}

absl::StatusOr<std::vector<std::string>> Tokenize(
    absl::string_view cif_string) {
  std::vector<std::string> tokens;

  absl::Status tokenize_status =
      TokenizeInternal(cif_string, [&](absl::string_view token) {
        tokens.emplace_back(token);
        return absl::OkStatus();
      });
  if (!tokenize_status.ok()) {
    return tokenize_status;
  }
  return tokens;
}

absl::StatusOr<std::vector<absl::string_view>> SplitLine(
    absl::string_view line) {
  std::vector<absl::string_view> tokens;
  absl::Status split_status =
      SplitLineInline(line, [&](absl::string_view token) {
        tokens.push_back(token);
        return absl::OkStatus();
      });
  if (!split_status.ok()) {
    return split_status;
  }
  return tokens;
}

absl::StatusOr<absl::flat_hash_map<std::string, CifDict>> ParseMultiDataCifDict(
    absl::string_view cif_string) {
  absl::flat_hash_map<std::string, CifDict> mapping;
  constexpr absl::string_view delimitor = "data_";
  // Check cif_string starts with correct offset.
  if (!cif_string.empty() && !absl::StartsWith(cif_string, delimitor)) {
    return absl::InvalidArgumentError(
        "Invalid format. MultiDataCifDict must start with 'data_'");
  }
  for (absl::string_view data_block :
       absl::StrSplit(cif_string, delimitor, absl::SkipEmpty())) {
    absl::string_view block_with_delimitor(
        data_block.data() - delimitor.size(),
        data_block.size() + delimitor.size());
    absl::StatusOr<CifDict> parsed_block =
        CifDict::FromString(block_with_delimitor);
    if (!parsed_block.ok()) {
      return parsed_block.status();
    }
    absl::string_view data_name = parsed_block->GetDataName();
    mapping[data_name] = *std::move(parsed_block);
  }

  return mapping;
}

}  // namespace alphafold3
