// Copyright 2024 DeepMind Technologies Limited
//
// AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
// this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// To request access to the AlphaFold 3 model parameters, follow the process set
// out at https://github.com/google-deepmind/alphafold3. You may only use these
// if received directly from Google. Use is subject to terms of use available at
// https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

#include <Python.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <new>
#include <optional>
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

class ProcessResidue {
 public:
  explicit ProcessResidue(const char* residue)
      : residue_(PyUnicode_InternFromString(residue)) {}
  bool IsResidue(PyObject* residue) {
    return ArePyObjectsEqual(residue_.get(), residue);
  }

  static bool ArePyObjectsEqual(PyObject* lhs, PyObject* rhs) {
    switch (PyObject_RichCompareBool(lhs, rhs, Py_EQ)) {
      case -1:
        PyErr_Clear();
        return false;
      case 0:
        return false;
      default:
        return true;
    }
  }

 private:
  ScopedPyObject residue_;
};

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

class FixArginine : public ProcessResidue {
 public:
  FixArginine()
      : ProcessResidue("ARG"),
        cd_(PyUnicode_InternFromString("CD")),
        nh1_(PyUnicode_InternFromString("NH1")),
        nh2_(PyUnicode_InternFromString("NH2")),
        hh11_(PyUnicode_InternFromString("HH11")),
        hh21_(PyUnicode_InternFromString("HH21")),
        hh12_(PyUnicode_InternFromString("HH12")),
        hh22_(PyUnicode_InternFromString("HH22")) {}
  void Fix(absl::Span<PyObject*> atom_ids, absl::Span<const float> atom_x,
           absl::Span<const float> atom_y, absl::Span<const float> atom_z) {
    std::ptrdiff_t cd_index = -1;
    std::ptrdiff_t nh1_index = -1;
    std::ptrdiff_t nh2_index = -1;
    std::ptrdiff_t hh11_index = -1;
    std::ptrdiff_t hh21_index = -1;
    std::ptrdiff_t hh12_index = -1;
    std::ptrdiff_t hh22_index = -1;
    for (std::ptrdiff_t index = 0; index < atom_ids.size(); ++index) {
      PyObject* atom_id = atom_ids[index];
      if (cd_index == -1 && ArePyObjectsEqual(atom_id, cd_.get())) {
        cd_index = index;
      } else if (nh1_index == -1 && ArePyObjectsEqual(atom_id, nh1_.get())) {
        nh1_index = index;
      } else if (nh2_index == -1 && ArePyObjectsEqual(atom_id, nh2_.get())) {
        nh2_index = index;
      } else if (hh11_index == -1 && ArePyObjectsEqual(atom_id, hh11_.get())) {
        hh11_index = index;
      } else if (hh21_index == -1 && ArePyObjectsEqual(atom_id, hh21_.get())) {
        hh21_index = index;
      } else if (hh12_index == -1 && ArePyObjectsEqual(atom_id, hh12_.get())) {
        hh12_index = index;
      } else if (hh22_index == -1 && ArePyObjectsEqual(atom_id, hh22_.get())) {
        hh22_index = index;
      }
    }
    if (cd_index < 0 || nh1_index < 0 || nh2_index < 0) {
      return;
    }
    Position3 cd_pos(atom_x[cd_index], atom_y[cd_index], atom_z[cd_index]);
    Position3 nh1_pos(atom_x[nh1_index], atom_y[nh1_index], atom_z[nh1_index]);
    Position3 nh2_pos(atom_x[nh2_index], atom_y[nh2_index], atom_z[nh2_index]);
    if (DistanceSquared(nh1_pos, cd_pos) <= DistanceSquared(nh2_pos, cd_pos)) {
      return;
    }
    std::swap(atom_ids[nh1_index], atom_ids[nh2_index]);
    if (hh11_index >= 0 && hh21_index >= 0) {
      std::swap(atom_ids[hh11_index], atom_ids[hh21_index]);
    } else if (hh11_index >= 0) {
      Py_DECREF(atom_ids[hh11_index]);
      Py_INCREF(hh21_.get());
      atom_ids[hh11_index] = hh21_.get();
    } else if (hh21_index >= 0) {
      Py_DECREF(atom_ids[hh21_index]);
      Py_INCREF(hh11_.get());
      atom_ids[hh21_index] = hh11_.get();
    }
    if (hh12_index >= 0 && hh22_index >= 0) {
      std::swap(atom_ids[hh12_index], atom_ids[hh22_index]);
    } else if (hh12_index >= 0) {
      Py_DECREF(atom_ids[hh12_index]);
      Py_INCREF(hh22_.get());
      atom_ids[hh12_index] = hh22_.get();
    } else if (hh22_index >= 0) {
      Py_DECREF(atom_ids[hh22_index]);
      Py_INCREF(hh21_.get());
      atom_ids[hh22_index] = hh21_.get();
    }
  }

 private:
  ScopedPyObject cd_;
  ScopedPyObject nh1_;
  ScopedPyObject nh2_;
  ScopedPyObject hh11_;
  ScopedPyObject hh21_;
  ScopedPyObject hh12_;
  ScopedPyObject hh22_;
};

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

std::pair<py::object, MmcifLayout> MmcifFilter(  //
    const CifDict& mmcif,                        //
    bool include_nucleotides,                    //
    bool include_ligands,                        //
    bool include_water,                          //
    bool include_other,                          //
    absl::string_view model_id) {
  if (_import_array() < 0) {
    throw py::import_error("Failed to import NumPy.");
  }
  auto layout = ReadMmcifLayout(mmcif, model_id);
  std::unique_ptr<std::vector<std::uint64_t>> keep_indices;
  size_t new_num_atoms;

  {
    py::gil_scoped_release release;

    AtomSiteLoop atom_site(mmcif);

    auto keep_chain_ids =
        SelectChains(mmcif, include_nucleotides, include_ligands, include_water,
                     include_other);

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
      std::vector<std::size_t> indices(keep_indices->begin(), keep_indices->end());
      std::size_t* start = indices.data();
      std::size_t num_atom = indices.size();
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

  return std::make_pair(py::reinterpret_steal<py::object>(arr),
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
  std::optional<FixArginine> arginine;
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

  if (fix_arginine) {
    arginine.emplace();
  }
  if (!arginine.has_value()) {
    return;
  }

  for (std::size_t res_index = 0; res_index < layout.num_residues();
       ++res_index) {
    auto [atom_start, atom_end] = layout.atom_range(res_index);
    std::size_t atom_count = atom_end - atom_start;
    PyObject* resname = comp_id[atom_start];
    if (arginine.has_value() && arginine->IsResidue(resname)) {
      arginine->Fix(atom_id.subspan(atom_start, atom_count),
                    atom_x.subspan(atom_start, atom_count),
                    atom_y.subspan(atom_start, atom_count),
                    atom_z.subspan(atom_start, atom_count));
    }
  }
}

std::vector<bool> SelectedPolymerResidueMask(
    const MmcifLayout& layout,
    const std::vector<absl::string_view>& atom_site_label_asym_ids,  //
    const std::vector<absl::string_view>& atom_site_label_seq_ids,   //
    const std::vector<absl::string_view>& atom_site_label_comp_ids,  //
    const std::vector<absl::string_view>& poly_seq_asym_ids,         //
    const std::vector<absl::string_view>& poly_seq_seq_ids,          //
    const std::vector<absl::string_view>& poly_seq_mon_ids           //
) {
  absl::flat_hash_map<std::pair<absl::string_view, absl::string_view>,
                      absl::string_view>
      selected;
  selected.reserve(layout.num_residues());
  // layout.residues() is O(1) while layout.residue_starts() is O(num_res).
  const std::vector<std::size_t>& residue_starts = layout.residue_starts();
  for (int i = 0; i < layout.residues().size(); ++i) {
    std::size_t res_start = residue_starts[i];
    std::size_t res_end = layout.residues()[i];
    if (res_start == res_end) {
      continue;  // Skip empty residues (containing no atoms).
    }

    absl::string_view label_seq_id = atom_site_label_seq_ids[i];
    if (label_seq_id == ".") {
      continue;  // Skip non-polymers.
    }

    absl::string_view label_asym_id = atom_site_label_asym_ids[i];
    absl::string_view label_comp_id = atom_site_label_comp_ids[i];
    selected[std::make_pair(label_asym_id, label_seq_id)] = label_comp_id;
  }

  std::vector<bool> mask;
  mask.reserve(poly_seq_mon_ids.size());
  for (int i = 0; i < poly_seq_mon_ids.size(); ++i) {
    absl::string_view poly_seq_asym_id = poly_seq_asym_ids[i];
    absl::string_view poly_seq_seq_id = poly_seq_seq_ids[i];
    absl::string_view poly_seq_mon_id = poly_seq_mon_ids[i];

    auto it = selected.find(std::make_pair(poly_seq_asym_id, poly_seq_seq_id));
    if (it != selected.end()) {
      mask.push_back(it->second == poly_seq_mon_id);
    } else {
      mask.push_back(true);  // Missing residues are not heterogeneous.
    }
  }
  return mask;
}

std::pair<std::vector<bool>, std::vector<bool>> SelectedLigandResidueMask(
    const MmcifLayout& layout,                                           //
    const std::vector<absl::string_view>& atom_site_label_asym_ids,      //
    const std::vector<absl::string_view>& atom_site_label_seq_ids,       //
    const std::vector<absl::string_view>& atom_site_auth_seq_ids,        //
    const std::vector<absl::string_view>& atom_site_label_comp_ids,      //
    const std::vector<absl::string_view>& atom_site_pdbx_pdb_ins_codes,  //
    const std::vector<absl::string_view>& nonpoly_asym_ids,              //
    const std::vector<absl::string_view>& nonpoly_auth_seq_ids,          //
    const std::vector<absl::string_view>& nonpoly_pdb_ins_codes,         //
    const std::vector<absl::string_view>& nonpoly_mon_ids,               //
    const std::vector<absl::string_view>& branch_asym_ids,               //
    const std::vector<absl::string_view>& branch_auth_seq_ids,           //
    const std::vector<absl::string_view>& branch_pdb_ins_codes,          //
    const std::vector<absl::string_view>& branch_mon_ids) {
  absl::flat_hash_map<
      std::tuple<absl::string_view, absl::string_view, absl::string_view>,
      absl::string_view>
      selected;
  selected.reserve(layout.num_residues());
  // layout.residues() is O(1) while layout.residue_starts() is O(num_res).
  const std::vector<std::size_t>& residue_starts = layout.residue_starts();
  for (int i = 0; i < layout.residues().size(); ++i) {
    std::size_t res_start = residue_starts[i];
    std::size_t res_end = layout.residues()[i];
    if (res_start == res_end) {
      continue;  // Skip empty residues (containing no atoms).
    }

    absl::string_view label_seq_id = atom_site_label_seq_ids[i];
    if (label_seq_id != ".") {
      continue;  // Skip polymers.
    }

    absl::string_view label_asym_id = atom_site_label_asym_ids[i];
    absl::string_view auth_seq_id = atom_site_auth_seq_ids[i];
    absl::string_view ins_code = atom_site_pdbx_pdb_ins_codes[i];
    ins_code = ins_code == "?" ? "." : ins_code;  // Remap unknown to unset.
    absl::string_view label_comp_id = atom_site_label_comp_ids[i];
    selected[std::make_tuple(label_asym_id, auth_seq_id, ins_code)] =
        label_comp_id;
  }

  std::vector<bool> nonpoly_mask;
  nonpoly_mask.reserve(nonpoly_asym_ids.size());
  for (int i = 0; i < nonpoly_asym_ids.size(); ++i) {
    absl::string_view nonpoly_asym_id = nonpoly_asym_ids[i];
    absl::string_view nonpoly_auth_seq_id = nonpoly_auth_seq_ids[i];
    absl::string_view nonpoly_ins_code = nonpoly_pdb_ins_codes[i];
    // Remap unknown to unset.
    nonpoly_ins_code = nonpoly_ins_code == "?" ? "." : nonpoly_ins_code;
    absl::string_view nonpoly_mon_id = nonpoly_mon_ids[i];

    auto it = selected.find(std::make_tuple(
        nonpoly_asym_id, nonpoly_auth_seq_id, nonpoly_ins_code));
    if (it != selected.end()) {
      nonpoly_mask.push_back(it->second == nonpoly_mon_id);
    } else {
      nonpoly_mask.push_back(true);  // Missing residues are not heterogeneous.
    }
  }

  std::vector<bool> branch_mask;
  branch_mask.reserve(branch_asym_ids.size());
  for (int i = 0; i < branch_asym_ids.size(); ++i) {
    absl::string_view branch_asym_id = branch_asym_ids[i];
    absl::string_view branch_auth_seq_id = branch_auth_seq_ids[i];

    // Insertion codes in _pdbx_branch_scheme are not required and can be
    // missing. Default to unset ('.') in such case.
    absl::string_view branch_ins_code;
    if (i < branch_pdb_ins_codes.size()) {
      branch_ins_code = branch_pdb_ins_codes[i];
      // Remap unknown to unset.
      branch_ins_code = branch_ins_code == "?" ? "." : branch_ins_code;
    } else {
      branch_ins_code = ".";
    }

    absl::string_view branch_mon_id = branch_mon_ids[i];

    auto it = selected.find(
        std::make_tuple(branch_asym_id, branch_auth_seq_id, branch_ins_code));
    if (it != selected.end()) {
      branch_mask.push_back(it->second == branch_mon_id);
    } else {
      branch_mask.push_back(true);  // Missing residues are not heterogeneous.
    }
  }

  return std::make_pair(nonpoly_mask, branch_mask);
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

  m.def("selected_polymer_residue_mask", SelectedPolymerResidueMask,
        py::arg("layout"),                         //
        py::arg("atom_site_label_asym_ids"),       //
        py::arg("atom_site_label_seq_ids"),        //
        py::arg("atom_site_label_comp_ids"),       //
        py::arg("poly_seq_asym_ids"),              //
        py::arg("poly_seq_seq_ids"),               //
        py::arg("poly_seq_mon_ids"),               //
        py::call_guard<py::gil_scoped_release>(),  //
        py::doc(kSelectedPolymerResidueMask + 1)   //
  );

  m.def("selected_ligand_residue_mask", SelectedLigandResidueMask,
        py::arg("layout"),                         //
        py::arg("atom_site_label_asym_ids"),       //
        py::arg("atom_site_label_seq_ids"),        //
        py::arg("atom_site_auth_seq_ids"),         //
        py::arg("atom_site_label_comp_ids"),       //
        py::arg("atom_site_pdbx_pdb_ins_codes"),   //
        py::arg("nonpoly_asym_ids"),               //
        py::arg("nonpoly_auth_seq_ids"),           //
        py::arg("nonpoly_pdb_ins_codes"),          //
        py::arg("nonpoly_mon_ids"),                //
        py::arg("branch_asym_ids"),                //
        py::arg("branch_auth_seq_ids"),            //
        py::arg("branch_pdb_ins_codes"),           //
        py::arg("branch_mon_ids"),                 //
        py::call_guard<py::gil_scoped_release>(),  //
        py::doc(kSelectedLigandResidueMask + 1)    //
  );
}

}  // namespace alphafold3
