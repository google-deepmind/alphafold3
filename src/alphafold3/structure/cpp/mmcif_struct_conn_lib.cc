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

#include <cstddef>
#include <iterator>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "absl/algorithm/container.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "alphafold3/parsers/cpp/cif_dict_lib.h"
#include "alphafold3/structure/cpp/mmcif_struct_conn.h"

namespace alphafold3 {

namespace {

struct AtomId {
  absl::string_view chain_id;
  absl::string_view res_id_1;
  absl::string_view res_id_2;
  absl::string_view atom_name;
  absl::string_view alt_id;

  friend bool operator==(const AtomId&, const AtomId&) = default;
  template <typename H>
  friend H AbslHashValue(H h, const AtomId& m) {
    return H::combine(std::move(h), m.chain_id, m.res_id_1, m.res_id_2,
                      m.atom_name, m.alt_id);
  }
};

using StringArrayRef = absl::Span<const std::string>;
using BondIndexByAtom = absl::flat_hash_map<AtomId, std::vector<std::size_t>>;
using BondAtomIndices = std::vector<std::size_t>;

// Returns whether each container is the same size.
template <typename C, typename... Cs>
bool AreSameSize(const C& c, const Cs&... cs) {
  return ((c.size() == cs.size()) && ...);
}

struct ColumnSpec {
  absl::string_view chain_id_col;
  absl::string_view res_id_1_col;
  absl::string_view res_id_2_col;
  absl::string_view atom_name_col;
  std::optional<absl::string_view> alt_id_col;  // Not used by OpenMM.
};

class AtomColumns {
 public:
  static absl::StatusOr<AtomColumns> Create(const CifDict& mmcif,
                                            const ColumnSpec& column_spec) {
    StringArrayRef chain_id = mmcif[column_spec.chain_id_col];
    StringArrayRef res_id_1 = mmcif[column_spec.res_id_1_col];
    StringArrayRef res_id_2 = mmcif[column_spec.res_id_2_col];
    StringArrayRef atom_name = mmcif[column_spec.atom_name_col];
    if (!AreSameSize(chain_id, res_id_1, res_id_2, atom_name)) {
      return absl::InvalidArgumentError(absl::StrCat(
          "Atom columns are not the same size. ",                       //
          "len(", column_spec.chain_id_col, ")=", chain_id.size(),      //
          ", len(", column_spec.res_id_1_col, ")=", res_id_1.size(),    //
          ", len(", column_spec.res_id_2_col, ")=", res_id_2.size(),    //
          ", len(", column_spec.atom_name_col, ")=", atom_name.size(),  //
          "."));
    }
    if (column_spec.alt_id_col.has_value()) {
      StringArrayRef alt_id = mmcif[*column_spec.alt_id_col];
      if (!AreSameSize(alt_id, chain_id)) {
        return absl::InvalidArgumentError(absl::StrCat(
            "Atom columns are not the same size. ",                   //
            "len(", column_spec.chain_id_col, ")=", chain_id.size(),  //
            ", len(", *column_spec.alt_id_col, ")=", alt_id.size(),   //
            "."));
      }
      return AtomColumns(chain_id, res_id_1, res_id_2, atom_name, alt_id,
                         column_spec);
    } else {
      return AtomColumns(chain_id, res_id_1, res_id_2, atom_name, std::nullopt,
                         column_spec);
    }
  }

  inline std::size_t size() const { return size_; }

  absl::string_view GetNormalizedAltId(const std::size_t index) const {
    constexpr absl::string_view kFullStop = ".";
    if (alt_id_.has_value()) {
      absl::string_view alt_id = (*alt_id_)[index];
      return alt_id == "?" ? kFullStop : alt_id;
    } else {
      return kFullStop;
    }
  }

  AtomId GetAtom(const std::size_t index) const {
    return {.chain_id = chain_id_[index],
            .res_id_1 = res_id_1_[index],
            .res_id_2 = res_id_2_[index],
            .atom_name = atom_name_[index],
            .alt_id = GetNormalizedAltId(index)};
  }

  std::string GetAtomString(const std::size_t index) const {
    std::string alt_id_col;
    if (column_spec_.alt_id_col.has_value()) {
      alt_id_col = *column_spec_.alt_id_col;
    } else {
      alt_id_col = "default label_alt_id";
    }
    return absl::StrCat(
        column_spec_.chain_id_col, "=", chain_id_[index], ", ",    //
        column_spec_.res_id_1_col, "=", res_id_1_[index], ", ",    //
        column_spec_.res_id_2_col, "=", res_id_2_[index], ", ",    //
        column_spec_.atom_name_col, "=", atom_name_[index], ", ",  //
        alt_id_col, "=", GetNormalizedAltId(index));               //
  }

 private:
  AtomColumns(StringArrayRef chain_id, StringArrayRef res_id_1,
              StringArrayRef res_id_2, StringArrayRef atom_name,
              std::optional<StringArrayRef> alt_id,
              const ColumnSpec& column_spec)
      : chain_id_(chain_id),
        res_id_1_(res_id_1),
        res_id_2_(res_id_2),
        atom_name_(atom_name),
        alt_id_(alt_id),
        column_spec_(column_spec),
        size_(chain_id.size()) {}
  StringArrayRef chain_id_;
  StringArrayRef res_id_1_;
  StringArrayRef res_id_2_;
  StringArrayRef atom_name_;
  std::optional<StringArrayRef> alt_id_;
  ColumnSpec column_spec_;
  std::size_t size_;
};

// Adds the atom index to any rows in the bond table involving that atom.
absl::Status FillInBondsForAtom(const BondIndexByAtom& bond_index_by_atom,
                                const AtomId& atom,
                                const std::size_t atom_index,
                                BondAtomIndices& bond_atom_indices) {
  if (auto bond_index_it = bond_index_by_atom.find(atom);
      bond_index_it != bond_index_by_atom.end()) {
    for (std::size_t bond_index : bond_index_it->second) {
      if (bond_index < 0 || bond_index >= bond_atom_indices.size()) {
        return absl::OutOfRangeError(
            absl::StrCat("Bond index out of range: ", bond_index));
      }
      bond_atom_indices[bond_index] = atom_index;
    }
  }
  return absl::OkStatus();
}

// Checks that the CifDict has all of the columns in the column spec.
bool HasAllColumns(const CifDict& mmcif, const ColumnSpec& columns) {
  return mmcif.Contains(columns.chain_id_col) &&
         mmcif.Contains(columns.res_id_1_col) &&
         mmcif.Contains(columns.res_id_2_col) &&
         mmcif.Contains(columns.atom_name_col) &&
         (!columns.alt_id_col.has_value() ||
          mmcif.Contains(*columns.alt_id_col));
}

// Fully specified ptnr1 atom.
constexpr ColumnSpec kStructConnPtnr1ColumnsFull{
    .chain_id_col = "_struct_conn.ptnr1_label_asym_id",
    .res_id_1_col = "_struct_conn.ptnr1_auth_seq_id",
    .res_id_2_col = "_struct_conn.pdbx_ptnr1_PDB_ins_code",
    .atom_name_col = "_struct_conn.ptnr1_label_atom_id",
    .alt_id_col = "_struct_conn.pdbx_ptnr1_label_alt_id",
};

// Fully specified ptnr2 atom.
constexpr ColumnSpec kStructConnPtnr2ColumnsFull{
    .chain_id_col = "_struct_conn.ptnr2_label_asym_id",
    .res_id_1_col = "_struct_conn.ptnr2_auth_seq_id",
    .res_id_2_col = "_struct_conn.pdbx_ptnr2_PDB_ins_code",
    .atom_name_col = "_struct_conn.ptnr2_label_atom_id",
    .alt_id_col = "_struct_conn.pdbx_ptnr2_label_alt_id",
};

// Columns used by OpenMM for ptnr1 atoms.
constexpr ColumnSpec kStructConnPtnr1OpenMM{
    .chain_id_col = "_struct_conn.ptnr1_label_asym_id",
    .res_id_1_col = "_struct_conn.ptnr1_label_seq_id",
    .res_id_2_col = "_struct_conn.ptnr1_label_comp_id",
    .atom_name_col = "_struct_conn.ptnr1_label_atom_id",
    .alt_id_col = std::nullopt,
};

// Columns used by OpenMM for ptnr2 atoms.
constexpr ColumnSpec kStructConnPtnr2OpenMM{
    .chain_id_col = "_struct_conn.ptnr2_label_asym_id",
    .res_id_1_col = "_struct_conn.ptnr2_label_seq_id",
    .res_id_2_col = "_struct_conn.ptnr2_label_comp_id",
    .atom_name_col = "_struct_conn.ptnr2_label_atom_id",
    .alt_id_col = std::nullopt,
};

// Fully specified atom sites.
constexpr ColumnSpec kAtomSiteColumnsFull{
    .chain_id_col = "_atom_site.label_asym_id",
    .res_id_1_col = "_atom_site.auth_seq_id",
    .res_id_2_col = "_atom_site.pdbx_PDB_ins_code",
    .atom_name_col = "_atom_site.label_atom_id",
    .alt_id_col = "_atom_site.label_alt_id",
};

// Atom site columns used to match OpenMM _struct_conn tables.
constexpr ColumnSpec kAtomSiteColumnsOpenMM{
    .chain_id_col = "_atom_site.label_asym_id",
    .res_id_1_col = "_atom_site.label_seq_id",
    .res_id_2_col = "_atom_site.label_comp_id",
    .atom_name_col = "_atom_site.label_atom_id",
    .alt_id_col = "_atom_site.label_alt_id",
};

}  // namespace

absl::StatusOr<std::pair<BondAtomIndices, BondAtomIndices>> GetBondAtomIndices(
    const CifDict& mmcif, absl::string_view model_id) {
  ColumnSpec ptnr1_columns, ptnr2_columns, atom_site_columns;

  if (HasAllColumns(mmcif, kStructConnPtnr1ColumnsFull) &&
      HasAllColumns(mmcif, kStructConnPtnr2ColumnsFull)) {
    ptnr1_columns = kStructConnPtnr1ColumnsFull;
    ptnr2_columns = kStructConnPtnr2ColumnsFull;
    atom_site_columns = kAtomSiteColumnsFull;
  } else {
    ptnr1_columns = kStructConnPtnr1OpenMM;
    ptnr2_columns = kStructConnPtnr2OpenMM;
    atom_site_columns = kAtomSiteColumnsOpenMM;
  }

  absl::StatusOr<AtomColumns> ptnr1_atoms =
      AtomColumns::Create(mmcif, ptnr1_columns);
  if (!ptnr1_atoms.ok()) {
    return ptnr1_atoms.status();
  }
  absl::StatusOr<AtomColumns> ptnr2_atoms =
      AtomColumns::Create(mmcif, ptnr2_columns);
  if (!ptnr2_atoms.ok()) {
    return ptnr2_atoms.status();
  }
  StringArrayRef struct_conn_id = mmcif["_struct_conn.id"];
  if (!AreSameSize(struct_conn_id, *ptnr1_atoms, *ptnr2_atoms)) {
    return absl::InvalidArgumentError(absl::StrCat(
        "Invalid '_struct_conn.' loop. ",                  //
        "len(id) = ", struct_conn_id.size(), ", ",         //
        "len(ptnr1_atoms) = ", ptnr1_atoms->size(), ", ",  //
        "len(ptnr2_atoms) = ", ptnr2_atoms->size(), "."    //
        ));
  }

  absl::StatusOr<AtomColumns> atoms =
      AtomColumns::Create(mmcif, atom_site_columns);
  if (!atoms.ok()) {
    return atoms.status();
  }
  StringArrayRef atom_site_id = mmcif["_atom_site.id"];
  StringArrayRef atom_site_model_id = mmcif["_atom_site.pdbx_PDB_model_num"];
  if (!AreSameSize(atom_site_id, atom_site_model_id, *atoms)) {
    return absl::InvalidArgumentError(absl::StrCat(
        "Invalid '_atom_site.' loop. ",                                //
        "len(id)= ", atom_site_id.size(), ", ",                        //
        "len(pdbx_PDB_model_num)= ", atom_site_model_id.size(), ", ",  //
        "len(atoms)= ", atoms->size(), "."));                          //
  }

  // Build maps from atom ID tuples to the rows in _struct_conn where that
  // atom appears (NB could be multiple).
  const std::size_t struct_conn_size = struct_conn_id.size();
  BondIndexByAtom ptnr1_rows_by_atom(struct_conn_size);
  BondIndexByAtom ptnr2_rows_by_atom(struct_conn_size);
  for (std::size_t i = 0; i < struct_conn_size; ++i) {
    ptnr1_rows_by_atom[ptnr1_atoms->GetAtom(i)].push_back(i);
    ptnr2_rows_by_atom[ptnr2_atoms->GetAtom(i)].push_back(i);
  }

  // Allocate two output arrays with one element per row in struct_conn, where
  // each element will be the index of that atom in the atom_site table.
  // Fill the arrays with atom_site_size, which is an invalid value, so that
  // we can check at the end that each atom has been found.
  const std::size_t atom_site_size = atom_site_id.size();
  BondAtomIndices ptnr1_atom_indices(struct_conn_size, atom_site_size);
  BondAtomIndices ptnr2_atom_indices(struct_conn_size, atom_site_size);

  bool model_id_ecountered = false;
  absl::flat_hash_set<absl::string_view> seen_alt_ids;
  for (std::size_t atom_i = 0; atom_i < atom_site_size; ++atom_i) {
    if (atom_site_model_id[atom_i] != model_id) {
      if (!model_id_ecountered) {
        continue;
      } else {
        // Models are contiguous so once we see a different model ID after
        // encountering our model ID then we can exit early.
        break;
      }
    } else {
      model_id_ecountered = true;
    }
    AtomId atom = atoms->GetAtom(atom_i);
    seen_alt_ids.insert(atom.alt_id);

    if (auto fill_in_bonds_status1 = FillInBondsForAtom(
            ptnr1_rows_by_atom, atom, atom_i, ptnr1_atom_indices);
        !fill_in_bonds_status1.ok()) {
      return fill_in_bonds_status1;
    }
    if (auto fill_in_bonds_status2 = FillInBondsForAtom(
            ptnr2_rows_by_atom, atom, atom_i, ptnr2_atom_indices);
        !fill_in_bonds_status2.ok()) {
      return fill_in_bonds_status2;
    }
  }
  // The seen_alt_ids check is a workaround for a known PDB issue: some mmCIFs
  // (2evw, 2g0v, 2g0x, 2g0z, 2g10, 2g11, 2g12, 2g14, 2grz, 2ntw as of 2024)
  // have multiple models and they set different whole-chain altloc in each
  // model. The bond table however doesn't distinguish between models, so there
  // are bonds that are valid only for some models. E.g. 2grz has model 1 with
  // chain A with altloc A, and model 2 with chain A with altloc B. The bonds
  // table lists a bond for each of these.

  // Check that a ptnr1 atom was found for every bond.
  if (auto row_it = absl::c_find(ptnr1_atom_indices, atom_site_size);
      row_it != ptnr1_atom_indices.end()) {
    if (seen_alt_ids.size() > 1 || seen_alt_ids.contains(".") ||
        seen_alt_ids.contains("?")) {
      std::size_t i = std::distance(ptnr1_atom_indices.begin(), row_it);
      return absl::InvalidArgumentError(
          absl::StrCat("Error parsing \"", mmcif.GetDataName(), "\". ",
                       "Cannot find atom for bond ID ", struct_conn_id[i], ": ",
                       ptnr1_atoms->GetAtomString(i)));
    }
  }

  // Check that a ptnr2 atom was found for every bond.
  if (auto row_it = absl::c_find(ptnr2_atom_indices, atom_site_size);
      row_it != ptnr2_atom_indices.end()) {
    if (seen_alt_ids.size() > 1 || seen_alt_ids.contains(".") ||
        seen_alt_ids.contains("?")) {
      std::size_t i = std::distance(ptnr2_atom_indices.begin(), row_it);
      return absl::InvalidArgumentError(
          absl::StrCat("Error parsing \"", mmcif.GetDataName(), "\". ",
                       "Cannot find atom for bond ID ", struct_conn_id[i], ": ",
                       ptnr2_atoms->GetAtomString(i)));
    }
  }

  if (!model_id_ecountered) {
    return absl::InvalidArgumentError(absl::StrCat(
        "Error parsing \"", mmcif.GetDataName(), "\". model_id \"", model_id,
        "\" not found in _atom_site.pdbx_PDB_model_num."));
  }

  return std::make_pair(std::move(ptnr1_atom_indices),
                        std::move(ptnr2_atom_indices));
}

}  // namespace alphafold3
