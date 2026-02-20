// Copyright 2024 DeepMind Technologies Limited
//
// AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
// this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// To request access to the AlphaFold 3 model parameters, follow the process set
// out at https://github.com/google-deepmind/alphafold3. You may only use these
// if received directly from Google. Use is subject to terms of use available at
// https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

#include "alphafold3/model/mkdssp_pybind.h"

#include <filesystem>

#include <cif++/file.hpp>
#include <cif++/pdb.hpp>
#include <dssp.hpp>
#include <sstream>

#include "absl/strings/string_view.h"
#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"

namespace alphafold3 {
namespace py = pybind11;

void RegisterModuleMkdssp(pybind11::module m) {
  py::module site = py::module::import("site");
  py::list paths = py::cast<py::list>(site.attr("getsitepackages")());
  // Find the first path that contains the libcifpp components.cif file.
  bool found = false;
  for (const auto& py_path : paths) {
    auto path_str =
        std::filesystem::path(py::cast<absl::string_view>(py_path)) /
        "share/libcifpp/components.cif";
    if (std::filesystem::exists(path_str)) {
      setenv("LIBCIFPP_DATA_DIR", path_str.parent_path().c_str(), 0);
      found = true;
      break;
    }
  }
  if (!found) {
    // Fallback: check sysconfig data path for non-prefix installations.
    py::module sysconfig = py::module::import("sysconfig");
    auto data_path =
        std::filesystem::path(
            py::cast<absl::string_view>(sysconfig.attr("get_path")("data"))) /
        "share/libcifpp/components.cif";
    if (std::filesystem::exists(data_path)) {
      setenv("LIBCIFPP_DATA_DIR", data_path.parent_path().c_str(), 0);
      found = true;
    }
  }
  if (!found) {
    throw py::type_error("Could not find the libcifpp components.cif file.");
  }
  m.def(
      "get_dssp",
      [](absl::string_view mmcif, int model_no,
         int min_poly_proline_stretch_length,
         bool calculate_surface_accessibility) {
        cif::file cif_file(mmcif.data(), mmcif.size());
        dssp result(cif_file.front(), model_no, min_poly_proline_stretch_length,
                    calculate_surface_accessibility);
        std::stringstream sstream;
        result.write_legacy_output(sstream);
        return sstream.str();
      },
      py::arg("mmcif"), py::arg("model_no") = 1,
      py::arg("min_poly_proline_stretch_length") = 3,
      py::arg("calculate_surface_accessibility") = false,
      py::doc("Gets secondary structure from an mmCIF file."));
}

}  // namespace alphafold3
