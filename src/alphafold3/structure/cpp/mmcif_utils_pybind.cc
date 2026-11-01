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

#include <Python.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <new>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "numpy/ndarrayobject.h"
#include "numpy/ndarraytypes.h"
#include "numpy/npy_common.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/memory/memory.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "alphafold3/parsers/cpp/cif_dict_lib.h"
#include "alphafold3/structure/cpp/mmcif_altlocs.h"
#include "alphafold3/structure/cpp/mmcif_layout.h"
#include "pybind11/cast.h"
#include "pybind11/gil.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"
#include "pybind11/stl.h"
#include "pybind11_abseil/absl_casters.h"

namespace alphafold3 {
namespace {
namespace py = pybind11;

struct PyObjectDeleter {
  inline void operator()(PyObject* obj) const { Py_CLEAR(obj); }
};

using ScopedPyObject = std::unique_ptr<PyObject, PyObjectDeleter>;

using StringArrayRef = absl::Span<const std::string>;
using Indexer = absl::flat_hash_map<absl::string_view, std::size_t>;

// Returns the reverse look-up map of name to index.
Indexer MakeIndex(StringArrayRef col) {
  Indexer index;
  index.reserve(col.size());
  for (std::size_t i = 0; i < col.size(); ++i) {
    index[col[i]] = i;
  }
  return index;
}

// Returns whether each container is the same size.
template <typename C, typename... Cs>
bool AreSameSize(C c, const Cs&... cs) {
  return ((c.size() == cs.size()) && ...);
}

// Stores references to columns in `_atom_site` ensuring they all exist and
// are the same size.
struct AtomSiteLoop {
  explicit AtomSiteLoop(const CifDict& cif_dict)
      : id(cif_dict["_atom_site.id"]),
        model_id(cif_dict["_atom_site.pdbx_PDB_model_num"]),
        chain_id(cif_dict["_atom_site.label_asym_id"]),
        seq_id(cif_dict["_atom_site.label_seq_id"]),

        comp_id(cif_dict["_atom_site.label_comp_id"]),
        atom_id(cif_dict["_atom_site.label_atom_id"]),

        alt_id(cif_dict["_atom_site.label_alt_id"]),
        occupancy(cif_dict["_atom_site.occupancy"])

  {
    if (!AreSameSize(id, model_id, chain_id, seq_id, comp_id, atom_id, alt_id,
                     occupancy)) {
      throw py::value_error(
          absl::StrCat("Invalid '_atom_site.' loop. ",                     //
                       "len(id)=", id.size(), ", ",                        //
                       "len(pdbx_PDB_model_num)=", model_id.size(), ", ",  //
                       "len(label_asym_id)=", chain_id.size(), ", ",       //
                       "len(label_seq_id)=", seq_id.size(), ", ",          //
                       "len(label_comp_id)=", comp_id.size(), ", ",        //
                       "len(atom_id)=", atom_id.size(), ", ",              //
                       "len(label_alt_id)=", alt_id.size(), ", ",          //
                       "len(occupancy)=", occupancy.size()));
    }
  }
  StringArrayRef id;
  StringArrayRef model_id;
  StringArrayRef chain_id;
  StringArrayRef seq_id;
  StringArrayRef comp_id;
  StringArrayRef atom_id;
  StringArrayRef alt_id;
  StringArrayRef occupancy;
};

// Stores references to columns in `_entity` ensuring they all exist and are the
// same size.
struct EntityLoop {
  explicit EntityLoop(const CifDict& cif_dict)
      : id(cif_dict["_entity.id"]), type(cif_dict["_entity.type"]) {
    if (!AreSameSize(id, type)) {
      throw py::value_error(absl::StrCat("Invalid '_entity.' loop. ",  //
                                         "len(id)=", id.size(), ", ",  //
                                         "len(type)=", type.size()));
    }
  }
  StringArrayRef id;
  StringArrayRef type;
};

// Stores references to columns in `_entity_poly` ensuring they all exist and
// are the same size.
struct EntityPolyLoop {
  explicit EntityPolyLoop(const CifDict& cif_dict)
      : entity_id(cif_dict["_entity_poly.entity_id"]),
        type(cif_dict["_entity_poly.type"]) {
    if (!AreSameSize(entity_id, type)) {
      throw py::value_error(absl::StrCat("Invalid '_entity_poly.' loop. ",  //
                                         "len(entity_id)=", entity_id.size(),
                                         ", ",  //
                                         "len(type)=", type.size()));
    }
  }
  StringArrayRef entity_id;
  StringArrayRef type;
};

// Returns a set of entity names removing ones not included by the flags
// specified.
absl::flat_hash_set<absl::string_view> SelectChains(const CifDict& mmcif,
                                                    bool include_nucleotides,
                                                    bool include_ligands,
                                                    bool include_water,
                                                    bool include_other) {
  EntityLoop entity_loop(mmcif);
  EntityPolyLoop entity_poly(mmcif);
  absl::flat_hash_set<absl::string_view> permitted_polymers{"polypeptide(L)"};
  absl::flat_hash_set<absl::string_view> forbidden_polymers;
  for (absl::string_view type :
       {"polydeoxyribonucleotide", "polyribonucleotide",
        "polydeoxyribonucleotide/polyribonucleotide hybrid"}) {
    if (include_nucleotides) {
      permitted_polymers.emplace(type);
    } else {
      forbidden_polymers.emplace(type);
    }
  }

  absl::flat_hash_set<absl::string_view> permitted_nonpoly_entity_types;
  absl::flat_hash_set<absl::string_view> forbidden_nonpoly_entity_types;
  for (absl::string_view type : {"non-polymer", "branched"}) {
    if (include_ligands) {
      permitted_nonpoly_entity_types.emplace(type);
    } else {
      forbidden_nonpoly_entity_types.emplace(type);
    }
  }
  absl::string_view water_type = "water";
  if (include_water) {
    permitted_nonpoly_entity_types.emplace(water_type);
  } else {
    forbidden_nonpoly_entity_types.emplace(water_type);
  }

  StringArrayRef chain_ids = mmcif["_struct_asym.id"];
  StringArrayRef entity_ids = mmcif["_struct_asym.entity_id"];
  Indexer chain_index = MakeIndex(chain_ids);
  Indexer entity_poly_index = MakeIndex(entity_poly.entity_id);
  Indexer entity_id_to_index = MakeIndex(entity_loop.id);

  absl::flat_hash_set<absl::string_view> keep_chain_id;
  for (std::size_t i = 0; i < chain_ids.size(); ++i) {
    absl::string_view chain_id = chain_ids[i];
    absl::string_view entity_id = entity_ids[i];
    if (entity_id_to_index.empty() ||
        entity_loop.type[entity_id_to_index[entity_id]] == "polymer") {
      if (auto it = entity_poly_index.find(entity_id);
          it != entity_poly_index.end()) {
        absl::string_view poly_type = entity_poly.type[it->second];
        if (include_other) {
          if (!forbidden_polymers.contains(poly_type)) {
            keep_chain_id.insert(chain_id);
          }
        } else {
          if (permitted_polymers.contains(poly_type)) {
            keep_chain_id.insert(chain_id);
          }
        }
      }
    } else {
      absl::string_view entity_type =
          entity_loop.type[entity_id_to_index[entity_id]];
      if (include_other) {
        if (!forbidden_nonpoly_entity_types.contains(entity_type)) {
          keep_chain_id.insert(chain_id);
          continue;
        }
      } else {
        if (permitted_nonpoly_entity_types.contains(entity_type)) {
          keep_chain_id.insert(chain_id);
          continue;
        }
      }
    }
  }
  return keep_chain_id;
}

struct Position3 {
  float x;
  float y;
  float z;
};

float DistanceSquared(Position3 v1, Position3 v2) {
  float dx = v1.x - v2.x;
  float dy = v1.y - v2.y;
  float dz = v1.z - v2.z;
  return dx * dx + dy * dy + dz * dz;
}

// Returns the layout of the mmCIF `_atom_site` table.
inline MmcifLayout ReadMmcifLayout(const CifDict& mmcif,
                                   absl::string_view model_id = "") {
  py::gil_scoped_release release;
  auto mmcif_layout = MmcifLayout::Create(mmcif, model_id);
  if (mmcif_layout.ok()) {
    return *mmcif_layout;
  }

  throw py::value_error(std::string(mmcif_layout.status().message()));
}
py::tuple MmcifFilter(         //
    const CifDict& mmcif,      //
    bool include_nucleotides,  //
    bool include_ligands,      //
    bool include_water,        //
    bool include_other,        //
    absl::string_view model_id) {
  if (_import_array() < 0) {
    throw py::import_error("Failed to import NumPy.");
  }
  auto layout = ReadMmcifLayout(mmcif, model_id);
  std::unique_ptr<std::vector<std::uint64_t>> keep_indices;
  size_t new_num_atoms;
  absl::flat_hash_set<absl::string_view> keep_chain_ids;

  {
    py::gil_scoped_release release;

    AtomSiteLoop atom_site(mmcif);

    keep_chain_ids = SelectChains(mmcif, include_nucleotides, include_ligands,
                                  include_water, include_other);

    std::vector<std::size_t> chain_indices;
    chain_indices.reserve(keep_chain_ids.size());
    for (std::size_t i = 0; i < layout.num_chains(); ++i) {
      if (keep_chain_ids.contains(
              atom_site.chain_id[layout.atom_site_from_chain_index(i)])) {
        chain_indices.push_back(i);
      }
    }

    keep_indices =
        absl::WrapUnique(new std::vector<std::uint64_t>(ResolveMmcifAltLocs(
            layout, atom_site.comp_id, atom_site.atom_id, atom_site.alt_id,
            atom_site.occupancy, chain_indices)));
    new_num_atoms = keep_indices->size();

    if (layout.num_models() > 1) {
      keep_indices->reserve(layout.num_models() * new_num_atoms);
      std::uint64_t* start = &(*keep_indices->begin());
      std::size_t num_atom = keep_indices->size();
      // Copy first model indices into all model indices offsetting each copy.
      for (std::size_t i = 1; i < layout.num_models(); ++i) {
        std::size_t offset = i * layout.num_atoms();
        std::transform(start, start + num_atom,
                       std::back_inserter(*keep_indices),
                       [offset](std::size_t v) { return v + offset; });
      }
    }
  }

  layout.Filter(*keep_indices);

  npy_intp shape[] = {static_cast<npy_intp>(layout.num_models()),
                      static_cast<npy_intp>(new_num_atoms)};
  PyObject* arr =
      PyArray_SimpleNewFromData(2, shape, NPY_INT64, keep_indices->data());
  // Create a capsule to hold the memory of the buffer so NumPy knows how to
  // delete it when done with it.
  PyObject* capsule = PyCapsule_New(
      keep_indices.release(), nullptr, +[](PyObject* capsule_cleanup) {
        void* memory = PyCapsule_GetPointer(capsule_cleanup, nullptr);
        delete static_cast<std::vector<std::size_t>*>(memory);
      });
  PyArray_SetBaseObject(reinterpret_cast<PyArrayObject*>(arr), capsule);

  return py::make_tuple(py::reinterpret_steal<py::object>(arr), keep_chain_ids,
                        std::move(layout));
}

void MmcifFixResidues(               //
    const MmcifLayout& layout,       //
    absl::Span<PyObject*> comp_id,   //
    absl::Span<PyObject*> atom_id,   //
    absl::Span<const float> atom_x,  //
    absl::Span<const float> atom_y,  //
    absl::Span<const float> atom_z,  //
    bool fix_arginine                //
) {
  std::size_t num_atoms = layout.num_atoms();
  if (comp_id.size() != num_atoms || atom_id.size() != num_atoms ||
      atom_x.size() != num_atoms || atom_y.size() != num_atoms ||
      atom_z.size() != num_atoms) {
    throw py::value_error(
        absl::StrCat("Sizes must match. ",                   //
                     "num_atoms=", num_atoms, ", ",          //
                     "len(comp_id)=", comp_id.size(), ", ",  //
                     "len(atom_id)=", atom_id.size(), ", ",  //
                     "len(atom_x)=", atom_x.size(), ", ",    //
                     "len(atom_y)=", atom_y.size(), ", ",    //
                     "len(atom_z)=", atom_z.size()));
  }

  if (!fix_arginine) {
    return;
  }

  auto is_ascii = [](PyObject* obj, const char* lit) {
    return PyUnicode_Check(obj) &&
           PyUnicode_CompareWithASCIIString(obj, lit) == 0;
  };
  ScopedPyObject hh11, hh21, hh12, hh22;
  auto replace_atom_id = [](PyObject*& slot, ScopedPyObject& cached,
                            const char* name) {
    if (!cached) {
      cached.reset(PyUnicode_InternFromString(name));
      if (!cached) {
        throw py::error_already_set();
      }
    }
    PyObject* old = std::exchange(slot, cached.get());
    Py_INCREF(slot);
    Py_DECREF(old);
  };

  for (std::size_t res_index = 0; res_index < layout.num_residues();
       ++res_index) {
    auto [atom_start, atom_end] = layout.atom_range(res_index);
    if (atom_start == atom_end || !is_ascii(comp_id[atom_start], "ARG")) {
      continue;
    }
    std::ptrdiff_t cd = -1, nh1 = -1, nh2 = -1;
    std::ptrdiff_t i_hh11 = -1, i_hh21 = -1, i_hh12 = -1, i_hh22 = -1;
    const std::ptrdiff_t atom_count = atom_end - atom_start;
    for (std::ptrdiff_t idx = 0; idx < atom_count; ++idx) {
      PyObject* a = atom_id[atom_start + idx];
      if (cd < 0 && is_ascii(a, "CD")) {
        cd = idx;
      } else if (nh1 < 0 && is_ascii(a, "NH1")) {
        nh1 = idx;
      } else if (nh2 < 0 && is_ascii(a, "NH2")) {
        nh2 = idx;
      } else if (i_hh11 < 0 && is_ascii(a, "HH11")) {
        i_hh11 = idx;
      } else if (i_hh21 < 0 && is_ascii(a, "HH21")) {
        i_hh21 = idx;
      } else if (i_hh12 < 0 && is_ascii(a, "HH12")) {
        i_hh12 = idx;
      } else if (i_hh22 < 0 && is_ascii(a, "HH22")) {
        i_hh22 = idx;
      }
    }
    if (cd < 0 || nh1 < 0 || nh2 < 0) {
      continue;
    }
    Position3 cd_pos{atom_x[atom_start + cd], atom_y[atom_start + cd],
                     atom_z[atom_start + cd]};
    Position3 nh1_pos{atom_x[atom_start + nh1], atom_y[atom_start + nh1],
                      atom_z[atom_start + nh1]};
    Position3 nh2_pos{atom_x[atom_start + nh2], atom_y[atom_start + nh2],
                      atom_z[atom_start + nh2]};
    if (DistanceSquared(nh1_pos, cd_pos) <= DistanceSquared(nh2_pos, cd_pos)) {
      continue;
    }
    PyObject** res_atoms = atom_id.data() + atom_start;
    std::swap(res_atoms[nh1], res_atoms[nh2]);
    if (i_hh11 >= 0 && i_hh21 >= 0) {
      std::swap(res_atoms[i_hh11], res_atoms[i_hh21]);
    } else if (i_hh11 >= 0) {
      replace_atom_id(res_atoms[i_hh11], hh21, "HH21");
    } else if (i_hh21 >= 0) {
      replace_atom_id(res_atoms[i_hh21], hh11, "HH11");
    }
    if (i_hh12 >= 0 && i_hh22 >= 0) {
      std::swap(res_atoms[i_hh12], res_atoms[i_hh22]);
    } else if (i_hh12 >= 0) {
      replace_atom_id(res_atoms[i_hh12], hh22, "HH22");
    } else if (i_hh22 >= 0) {
      replace_atom_id(res_atoms[i_hh22], hh12, "HH12");
    }
  }
}

using ObjArray =
    py::array_t<PyObject*, py::array::c_style | py::array::forcecast>;

// Returns views into the UTF-8 buffers of the `str` elements of `array`, valid
// for as long as `array` keeps its elements alive.
std::vector<absl::string_view> ExtractStringViews(const ObjArray& array,
                                                  std::size_t expected_size,
                                                  absl::string_view name) {
  if (array.size() != expected_size) {
    throw py::value_error(absl::StrCat(name, " must have length ",
                                       expected_size, ", got ", array.size()));
  }
  std::vector<absl::string_view> result;
  result.reserve(expected_size);
  for (PyObject* obj : absl::MakeConstSpan(array.data(), expected_size)) {
    Py_ssize_t size;
    const char* buf = PyUnicode_AsUTF8AndSize(obj, &size);
    if (buf == nullptr) throw py::error_already_set();
    result.emplace_back(buf, size);
  }
  return result;
}

py::array_t<int64_t> MapResidueKeys(
    const std::array<ObjArray, 4>& src_arrays,
    const std::array<ObjArray, 4>& query_arrays) {
  const std::size_t num_src = src_arrays[0].size();
  const std::size_t num_query = query_arrays[0].size();
  std::array<std::vector<absl::string_view>, 4> src, query;
  for (int c = 0; c < 4; ++c) {
    src[c] = ExtractStringViews(src_arrays[c], num_src, "src_arrays");
    query[c] = ExtractStringViews(query_arrays[c], num_query, "query_arrays");
  }
  py::array_t<int64_t> result(num_query);
  absl::Span<int64_t> out = absl::MakeSpan(result.mutable_data(), num_query);
  std::string missing_key;
  {
    py::gil_scoped_release release;
    using ResKey = std::array<absl::string_view, 4>;
    absl::flat_hash_map<ResKey, int64_t> index_by_key;
    index_by_key.reserve(num_src);
    for (std::size_t i = 0; i < num_src; ++i) {
      index_by_key[{src[0][i], src[1][i], src[2][i], src[3][i]}] = i;
    }
    for (std::size_t i = 0; i < num_query; ++i) {
      const ResKey key = {query[0][i], query[1][i], query[2][i], query[3][i]};
      auto it = index_by_key.find(key);
      if (it == index_by_key.end()) {
        missing_key = absl::StrCat("('", key[0], "', '", key[1], "', '",
                                   key[2], "', '", key[3], "')");
        break;
      }
      out[i] = it->second;
    }
  }
  if (!missing_key.empty()) throw py::key_error(missing_key);
  return result;
}

py::array_t<bool> SelectedPolymerResidueMask(
    const MmcifLayout& layout,
    const ObjArray& atom_site_label_asym_ids_arr,  //
    const ObjArray& atom_site_label_seq_ids_arr,   //
    const ObjArray& atom_site_label_comp_ids_arr,  //
    const ObjArray& poly_seq_asym_ids_arr,         //
    const ObjArray& poly_seq_seq_ids_arr,          //
    const ObjArray& poly_seq_mon_ids_arr           //
) {
  const std::size_t num_res = layout.residues().size();
  const std::size_t num_poly = poly_seq_asym_ids_arr.size();
  auto atom_site_label_asym_ids = ExtractStringViews(
      atom_site_label_asym_ids_arr, num_res, "atom_site_label_asym_ids");
  auto atom_site_label_seq_ids = ExtractStringViews(
      atom_site_label_seq_ids_arr, num_res, "atom_site_label_seq_ids");
  auto atom_site_label_comp_ids = ExtractStringViews(
      atom_site_label_comp_ids_arr, num_res, "atom_site_label_comp_ids");
  auto poly_seq_asym_ids =
      ExtractStringViews(poly_seq_asym_ids_arr, num_poly, "poly_seq_asym_ids");
  auto poly_seq_seq_ids =
      ExtractStringViews(poly_seq_seq_ids_arr, num_poly, "poly_seq_seq_ids");
  auto poly_seq_mon_ids =
      ExtractStringViews(poly_seq_mon_ids_arr, num_poly, "poly_seq_mon_ids");

  py::array_t<bool> mask(num_poly);
  bool* out = mask.mutable_data();
  {
    py::gil_scoped_release release;
    absl::flat_hash_map<std::pair<absl::string_view, absl::string_view>,
                        absl::string_view>
        selected;
    selected.reserve(layout.num_residues());
    // layout.residues() is O(1) while layout.residue_starts() is O(num_res).
    const std::vector<std::size_t>& residue_starts = layout.residue_starts();
    for (std::size_t i = 0; i < num_res; ++i) {
      if (residue_starts[i] == layout.residues()[i] ||
          atom_site_label_seq_ids[i] == ".") {
        continue;  // Skip empty residues and non-polymers.
      }
      selected[{atom_site_label_asym_ids[i], atom_site_label_seq_ids[i]}] =
          atom_site_label_comp_ids[i];
    }

    for (std::size_t i = 0; i < num_poly; ++i) {
      auto it = selected.find({poly_seq_asym_ids[i], poly_seq_seq_ids[i]});
      // Missing residues are not heterogeneous.
      out[i] = it == selected.end() || it->second == poly_seq_mon_ids[i];
    }
  }
  return mask;
}

std::pair<py::array_t<bool>, py::array_t<bool>> SelectedLigandResidueMask(
    const MmcifLayout& layout,                         //
    const ObjArray& atom_site_label_asym_ids_arr,      //
    const ObjArray& atom_site_label_seq_ids_arr,       //
    const ObjArray& atom_site_auth_seq_ids_arr,        //
    const ObjArray& atom_site_label_comp_ids_arr,      //
    const ObjArray& atom_site_pdbx_pdb_ins_codes_arr,  //
    const ObjArray& nonpoly_asym_ids_arr,              //
    const ObjArray& nonpoly_auth_seq_ids_arr,          //
    const ObjArray& nonpoly_pdb_ins_codes_arr,         //
    const ObjArray& nonpoly_mon_ids_arr,               //
    const ObjArray& branch_asym_ids_arr,               //
    const ObjArray& branch_auth_seq_ids_arr,           //
    const ObjArray& branch_pdb_ins_codes_arr,          //
    const ObjArray& branch_mon_ids_arr) {
  const std::size_t num_res = layout.residues().size();
  const std::size_t num_nonpoly = nonpoly_asym_ids_arr.size();
  const std::size_t num_branch = branch_asym_ids_arr.size();
  auto atom_site_label_asym_ids = ExtractStringViews(
      atom_site_label_asym_ids_arr, num_res, "atom_site_label_asym_ids");
  auto atom_site_label_seq_ids = ExtractStringViews(
      atom_site_label_seq_ids_arr, num_res, "atom_site_label_seq_ids");
  auto atom_site_auth_seq_ids = ExtractStringViews(
      atom_site_auth_seq_ids_arr, num_res, "atom_site_auth_seq_ids");
  auto atom_site_label_comp_ids = ExtractStringViews(
      atom_site_label_comp_ids_arr, num_res, "atom_site_label_comp_ids");
  auto atom_site_pdbx_pdb_ins_codes =
      ExtractStringViews(atom_site_pdbx_pdb_ins_codes_arr, num_res,
                         "atom_site_pdbx_pdb_ins_codes");
  auto nonpoly_asym_ids =
      ExtractStringViews(nonpoly_asym_ids_arr, num_nonpoly, "nonpoly_asym_ids");
  auto nonpoly_auth_seq_ids = ExtractStringViews(
      nonpoly_auth_seq_ids_arr, num_nonpoly, "nonpoly_auth_seq_ids");
  auto nonpoly_pdb_ins_codes = ExtractStringViews(
      nonpoly_pdb_ins_codes_arr, num_nonpoly, "nonpoly_pdb_ins_codes");
  auto nonpoly_mon_ids =
      ExtractStringViews(nonpoly_mon_ids_arr, num_nonpoly, "nonpoly_mon_ids");
  auto branch_asym_ids =
      ExtractStringViews(branch_asym_ids_arr, num_branch, "branch_asym_ids");
  auto branch_auth_seq_ids = ExtractStringViews(
      branch_auth_seq_ids_arr, num_branch, "branch_auth_seq_ids");
  // Insertion codes in _pdbx_branch_scheme are not required and can be
  // missing. Default to unset ('.') in such case.
  auto branch_pdb_ins_codes = ExtractStringViews(
      branch_pdb_ins_codes_arr,
      branch_pdb_ins_codes_arr.size() == 0 ? 0 : num_branch,
      "branch_pdb_ins_codes");
  auto branch_mon_ids =
      ExtractStringViews(branch_mon_ids_arr, num_branch, "branch_mon_ids");

  py::array_t<bool> nonpoly_mask(num_nonpoly);
  py::array_t<bool> branch_mask(num_branch);
  bool* nonpoly_out = nonpoly_mask.mutable_data();
  bool* branch_out = branch_mask.mutable_data();
  {
    py::gil_scoped_release release;
    absl::flat_hash_map<
        std::tuple<absl::string_view, absl::string_view, absl::string_view>,
        absl::string_view>
        selected;
    selected.reserve(layout.num_residues());
    // layout.residues() is O(1) while layout.residue_starts() is O(num_res).
    const std::vector<std::size_t>& residue_starts = layout.residue_starts();
    for (std::size_t i = 0; i < num_res; ++i) {
      if (residue_starts[i] == layout.residues()[i] ||
          atom_site_label_seq_ids[i] != ".") {
        continue;  // Skip empty residues and polymers.
      }
      absl::string_view ins_code = atom_site_pdbx_pdb_ins_codes[i];
      ins_code = ins_code == "?" ? "." : ins_code;  // Remap unknown to unset.
      selected[{atom_site_label_asym_ids[i], atom_site_auth_seq_ids[i],
                ins_code}] = atom_site_label_comp_ids[i];
    }

    for (std::size_t i = 0; i < num_nonpoly; ++i) {
      absl::string_view ins_code = nonpoly_pdb_ins_codes[i];
      ins_code = ins_code == "?" ? "." : ins_code;  // Remap unknown to unset.
      auto it = selected.find(
          {nonpoly_asym_ids[i], nonpoly_auth_seq_ids[i], ins_code});
      // Missing residues are not heterogeneous.
      nonpoly_out[i] = it == selected.end() || it->second == nonpoly_mon_ids[i];
    }

    for (std::size_t i = 0; i < num_branch; ++i) {
      absl::string_view ins_code =
          branch_pdb_ins_codes.empty() ? "." : branch_pdb_ins_codes[i];
      ins_code = ins_code == "?" ? "." : ins_code;  // Remap unknown to unset.
      auto it =
          selected.find({branch_asym_ids[i], branch_auth_seq_ids[i], ins_code});
      // Missing residues are not heterogeneous.
      branch_out[i] = it == selected.end() || it->second == branch_mon_ids[i];
    }
  }
  return {std::move(nonpoly_mask), std::move(branch_mask)};
}

constexpr char kReadMmcifLayout[] = R"(
Returns the layout of the cif_dict.

Args:
  mmcif: mmCIF to calculate the layout for.
  model_id: If non-empty the layout of the given model is returned
    otherwise the layout of all models are returned.
Raises:
  ValueError: if the mmCIF is malformed or the number of atoms in each
    model are inconsistent.
)";

constexpr char kMmcifFilter[] = R"(
Returns NumpyArray of selected rows in `_atom_site` and new layout.

Args:
  mmcif: mmCIF to filter.
  include_nucleotides: Whether to include polymer entities of type:
    "polypeptide(L)\", "polydeoxyribonucleotide", "polyribonucleotide".
    Otherwise only "polypeptide(L)\". ("polypeptide(D)\" is never included.)
  include_ligands: Whether to include non-polymer entities of type:
    "non-polymer", "branched".
  include_water: Whether to include entities of type water.
  include_other: Whether to include other (non-standard) entity types
    that are not covered by any of the above parameters.
  model_id: If non-empty the model with given name is selected otherwise
    all models are selected.

Returns:
  A tuple containing a numpy array with a shape (num_models, num_atoms)
  with the atom_site indices selected and the new layout.

Raises:
  ValueError error if mmCIF dict does not have all required fields.
)";

constexpr char kMmcifFixResidues[] = R"(
Fixes residue columns in-place.

Args:
  layout: layout from filter command.
  comp_id: '_atom_site.label_comp_id' of first model.
  group: '_atom_site.group_PDB' of first model.
  atom_id: '_atom_site.label_atom_id' of first model.
  type_symbol: '_atom_site.type_symbol' of first model.
  atom_x: '_atom_site.Cartn_x' of first model.
  atom_y: '_atom_site.Cartn_y' of first model.
  atom_z: '_atom_site.Cartn_z' of first model.
  fix_mse: Whether to convert MSE residues into MET residues.
  fix_arg: Whether to ensure the atoms in ARG are in the correct order.
  fix_unknown_dna: Whether to convert DNA residues from N to DN.
  dna_mask: Which atoms are from DNA chains.

Raises:
  ValueError: If shapes are invalid.
)";

constexpr char kSelectedPolymerResidueMask[] = R"(
Returns a _pdbx_poly_seq_scheme mask for selected hetero residues.

Should be called after filtering the layout using mmcif_utils.filter.

Args:
  layout: Layout defining the _atom_site residue selection.
  atom_site_label_asym_ids: Internal (label) chain ID, per selected residue.
  atom_site_label_seq_ids: Internal (label) residue ID, per selected residue.
  atom_site_label_comp_ids: Residue name, per selected residue.
  poly_seq_asym_ids: Internal (label) chain ID, per residue.
  poly_seq_seq_ids: Internal (label) residue ID, per residue.
  poly_seq_mon_ids: Residue name, per residue.

Returns:
  A mask for the _pdbx_poly_seq_scheme table. If residues are selected
  using this mask, they will have consistent heterogeneous residue
  selection with the _atom_site table.
)";

constexpr char kSelectedLigandResidueMask[] = R"(
Returns masks for selected ligand hetero residues.

Should be called after filtering the layout using mmcif_utils.filter.

Args:
  layout: Layout defining the _atom_site residue selection.
  atom_site_label_asym_ids: Internal (label) chain ID, per selected residue.
  atom_site_label_seq_ids: Internal (author) residue ID, per selected residue.
  atom_site_auth_seq_ids: External (author) residue ID, per selected residue.
  atom_site_label_comp_ids: Residue name, per selected residue.
  atom_site_pdbx_pdb_ins_codes: Insertion code, per selected residue.
  nonpoly_asym_ids: Internal (label) chain ID, per residue from
   _pdbx_nonpoly_scheme.
  nonpoly_auth_seq_ids: External (author) residue ID, per residue from
   _pdbx_nonpoly_scheme.
  nonpoly_pdb_ins_codes: Residue name, per residue from
   _pdbx_nonpoly_scheme.
  nonpoly_mon_ids: Insertion code, per residue from _pdbx_nonpoly_scheme.
  branch_asym_ids: Internal (label) chain ID, per residue from
   _pdbx_branch_scheme.
  branch_auth_seq_ids: External (author) residue ID, per residue from
   _pdbx_branch_scheme.
  branch_pdb_ins_codes: Residue name, per residue from _pdbx_branch_scheme.
  branch_mon_ids: Insertion code, per residue from _pdbx_branch_scheme.

Returns:
  A tuple with masks for _pdbx_nonpoly_scheme and _pdbx_branch_scheme. If
  residues are selected using these masks, they will have consistent
  heterogeneous residue selection with the _atom_site table.
)";

}  // namespace

void RegisterModuleMmcifUtils(pybind11::module m) {
  m.def("read_layout", ReadMmcifLayout,
        py::arg("mmcif"),              //
        py::arg("model_id") = "",      //
        py::doc(kReadMmcifLayout + 1)  //
  );

  m.def("filter", MmcifFilter,               //
        py::arg("mmcif"),                    //
        py::arg("include_nucleotides"),      //
        py::arg("include_ligands") = false,  //
        py::arg("include_water") = false,    //
        py::arg("include_other") = false,    //
        py::arg("model_id") = "",            //
        py::doc(kMmcifFilter + 1)            //
  );

  m.def("fix_residues", MmcifFixResidues,
        py::arg("layout"),              //
        py::arg("comp_id"),             //
        py::arg("atom_id"),             //
        py::arg("atom_x"),              //
        py::arg("atom_y"),              //
        py::arg("atom_z"),              //
        py::arg("fix_arg") = false,     //
        py::doc(kMmcifFixResidues + 1)  //
  );

  m.def("map_residue_keys", MapResidueKeys, py::arg("src_arrays"),
        py::arg("query_arrays"),
        "Maps 4-column residue keys in query_arrays to their row index in "
        "src_arrays.");

  m.def("selected_polymer_residue_mask", SelectedPolymerResidueMask,
        py::arg("layout"),                        //
        py::arg("atom_site_label_asym_ids"),      //
        py::arg("atom_site_label_seq_ids"),       //
        py::arg("atom_site_label_comp_ids"),      //
        py::arg("poly_seq_asym_ids"),             //
        py::arg("poly_seq_seq_ids"),              //
        py::arg("poly_seq_mon_ids"),              //
        py::doc(kSelectedPolymerResidueMask + 1)  //
  );

  m.def("selected_ligand_residue_mask", SelectedLigandResidueMask,
        py::arg("layout"),                            //
        py::arg("atom_site_label_asym_ids"),          //
        py::arg("atom_site_label_seq_ids"),           //
        py::arg("atom_site_auth_seq_ids"),            //
        py::arg("atom_site_label_comp_ids"),          //
        py::arg("atom_site_pdbx_pdb_ins_codes"),      //
        py::arg("nonpoly_asym_ids"),                  //
        py::arg("nonpoly_auth_seq_ids"),              //
        py::arg("nonpoly_pdb_ins_codes"),             //
        py::arg("nonpoly_mon_ids"),                   //
        py::arg("branch_asym_ids"),                   //
        py::arg("branch_auth_seq_ids"),               //
        py::arg("branch_pdb_ins_codes"),              //
        py::arg("branch_mon_ids"),                    //
        py::doc(kSelectedLigandResidueMask + 1)       //
  );
}

}  // namespace alphafold3
