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

#include <string>

#include "absl/strings/string_view.h"
#include "alphafold3/parsers/cpp/cif_dict_lib.h"
#include "alphafold3/structure/cpp/mmcif_struct_conn.h"
#include "pybind11/gil.h"
#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"
#include "pybind11/stl.h"

namespace alphafold3 {

namespace py = pybind11;

constexpr char kGetBondAtomIndices[] = R"(
Extracts the indices of the atoms that participate in bonds.

This function has a workaround for a known PDB issue: some mmCIFs have
(2evw, 2g0v, 2g0x, 2g0z, 2g10, 2g11, 2g12, 2g14, 2grz, 2ntw as of 2024)
multiple models and they set different whole-chain altloc in each model.
The bond table however doesn't distinguish between models, so there are
bonds that are valid only for some models. E.g. 2grz has model 1 with
chain A with altloc A, and model 2 with chain A with altloc B. The bonds
table lists a bond for each of these. This case is rather rare (10 cases
in PDB as of 2024). For the offending bonds, the returned atom index is
set to the size of the atom_site table, i.e. it is an invalid index.

Args:
  mmcif: The mmCIF object to process.
  model_id: The ID of the model that the returned atoms will belong to. This
    should be a value in the mmCIF's _atom_site.pdbx_PDB_model_num column.

Returns:
  Two lists of atom indices, `from_atoms` and `to_atoms`, each one having
  length num_bonds (as defined by _struct_conn, the bonds table). The bond
  i, defined by the i'th row in _struct_conn, is a bond from atom at index
  from_atoms[i], to the atom at index to_atoms[i]. The indices are simple
  0-based indexes into the columns of the _atom_site table in the input
  mmCIF, and do not necessarily correspond to the values in _atom_site.id,
  or any other column.
)";

void RegisterModuleMmcifStructConn(pybind11::module m) {
  m.def(
      "get_bond_atom_indices",
      [](const CifDict& mmcif, absl::string_view model_id) {
        auto result = GetBondAtomIndices(mmcif, model_id);
        if (result.ok()) {
          return *result;
        }
        throw py::value_error(std::string(result.status().message()));
      },
      py::arg("mmcif_dict"), py::arg("model_id"),
      py::doc(kGetBondAtomIndices + 1),
      py::call_guard<py::gil_scoped_release>());
}

}  // namespace alphafold3
