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

#include "absl/container/flat_hash_map.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "alphafold3/parsers/cpp/cif_dict_lib.h"
#include "pybind11/gil.h"
#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"
#include "pybind11/stl.h"
#include "pybind11_abseil/absl_casters.h"

namespace alphafold3 {
namespace {
namespace py = pybind11;

// If present, returns the _atom_site.type_symbol. If not, infers it using
// _atom_site.label_comp_id (residue name), _atom_site.label_atom_id (atom name)
// and the CCD.
py::list GetOrInferTypeSymbol(const CifDict& mmcif,
                              const py::object& atom_id_to_type_symbol) {
  const auto& type_symbol = mmcif["_atom_site.type_symbol"];
  const auto& id = mmcif["_atom_site.id"];
  if (!type_symbol.empty()) {
    if (!id.empty() && type_symbol.size() != id.size()) {
      throw py::value_error(
          absl::StrCat("Invalid '_atom_site.' loop. len(id)=", id.size(),
                       ", len(type_symbol)=", type_symbol.size()));
    }
    py::list patched_type_symbol(type_symbol.size());
    for (size_t i = 0; i < type_symbol.size(); ++i) {
      patched_type_symbol[i] = type_symbol[i];
    }
    return patched_type_symbol;
  }

  const auto& label_comp_id = mmcif["_atom_site.label_comp_id"];
  const auto& label_atom_id = mmcif["_atom_site.label_atom_id"];
  if (label_comp_id.size() != label_atom_id.size() ||
      (!id.empty() && id.size() != label_comp_id.size())) {
    throw py::value_error(
        absl::StrCat("Invalid '_atom_site.' loop. len(id)=", id.size(),
                     ", len(label_comp_id)=", label_comp_id.size(),
                     ", len(label_atom_id)=", label_atom_id.size()));
  }
  py::list patched_type_symbol(label_comp_id.size());
  absl::flat_hash_map<absl::string_view,
                      absl::flat_hash_map<absl::string_view, py::object>>
      cache;
  for (size_t i = 0; i < label_comp_id.size(); ++i) {
    auto& inner = cache[absl::string_view(label_comp_id[i])];
    auto [it, inserted] =
        inner.try_emplace(absl::string_view(label_atom_id[i]));
    if (inserted) {
      it->second = atom_id_to_type_symbol(label_comp_id[i], label_atom_id[i]);
    }
    patched_type_symbol[i] = it->second;
  }
  return patched_type_symbol;
}

absl::flat_hash_map<absl::string_view, absl::string_view>
GetInternalToAuthorChainIdMap(const CifDict& mmcif) {
  const auto& label_asym_ids = mmcif["_atom_site.label_asym_id"];
  const auto& auth_asym_ids = mmcif["_atom_site.auth_asym_id"];
  if (label_asym_ids.size() != auth_asym_ids.size()) {
    throw py::value_error(absl::StrCat(
        "Invalid '_atom_site.' loop. len(label_asym_id)=",
        label_asym_ids.size(), ", len(auth_asym_id)=", auth_asym_ids.size()));
  }

  absl::flat_hash_map<absl::string_view, absl::string_view> mapping;
  for (size_t i = 0, num_rows = label_asym_ids.size(); i < num_rows; ++i) {
    // Use only the first internal_chain_id occurrence to generate the mapping.
    // It should not matter as there should not be a case where a single
    // internal chain ID would map to more than one author chain IDs (i.e. the
    // mapping should be injective). Since we need this method to be fast, we
    // choose not to check it.
    mapping.emplace(label_asym_ids[i], auth_asym_ids[i]);
  }
  return mapping;
}

}  // namespace

namespace py = pybind11;

void RegisterModuleMmcifAtomSite(pybind11::module m) {
  m.def("get_or_infer_type_symbol", &GetOrInferTypeSymbol, py::arg("mmcif"),
        py::arg("atom_id_to_type_symbol"));

  m.def("get_internal_to_author_chain_id_map", &GetInternalToAuthorChainIdMap,
        py::arg("mmcif"), py::call_guard<py::gil_scoped_release>());
}

}  // namespace alphafold3
