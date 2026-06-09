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

#include "alphafold3/data/cpp/msa_profile_pybind.h"
#include "alphafold3/model/json_serialize_pybind.h"
#include "alphafold3/model/mkdssp_pybind.h"
#include "alphafold3/parsers/cpp/cif_dict_pybind.h"
#include "alphafold3/parsers/cpp/fasta_iterator_pybind.h"
#include "alphafold3/parsers/cpp/msa_conversion_pybind.h"
#include "alphafold3/structure/cpp/aggregation_pybind.h"
#include "alphafold3/structure/cpp/membership_pybind.h"
#include "alphafold3/structure/cpp/mmcif_atom_site_pybind.h"
#include "alphafold3/structure/cpp/mmcif_layout_pybind.h"
#include "alphafold3/structure/cpp/mmcif_struct_conn_pybind.h"
#include "alphafold3/structure/cpp/mmcif_utils_pybind.h"
#include "alphafold3/structure/cpp/string_array_pybind.h"
#include "pybind11/pybind11.h"

namespace alphafold3 {
namespace {

// Include all modules as submodules to simplify building.
PYBIND11_MODULE(cpp, m) {
  RegisterModuleCifDict(m.def_submodule("cif_dict"));
  RegisterModuleFastaIterator(m.def_submodule("fasta_iterator"));
  RegisterModuleMsaConversion(m.def_submodule("msa_conversion"));
  RegisterModuleMmcifLayout(m.def_submodule("mmcif_layout"));
  RegisterModuleMmcifStructConn(m.def_submodule("mmcif_struct_conn"));
  RegisterModuleMembership(m.def_submodule("membership"));
  RegisterModuleMmcifUtils(m.def_submodule("mmcif_utils"));
  RegisterModuleAggregation(m.def_submodule("aggregation"));
  RegisterModuleStringArray(m.def_submodule("string_array"));
  RegisterModuleMmcifAtomSite(m.def_submodule("mmcif_atom_site"));
  RegisterModuleJsonSerialize(m.def_submodule("json_serialize"));
  RegisterModuleMkdssp(m.def_submodule("mkdssp"));
  RegisterModuleMsaProfile(m.def_submodule("msa_profile"));
}

}  // namespace
}  // namespace alphafold3
