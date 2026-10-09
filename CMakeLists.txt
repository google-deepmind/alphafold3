# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

cmake_minimum_required(VERSION 3.28)

# This forces Git to use the 'files' backend for all FetchContent operations.
# This fixes libcifpp and dssp incompatibility with newer git versions.
set(ENV{GIT_CONFIG_PARAMETERS} "'init.defaultRefFormat=files'")

project(
  "${SKBUILD_PROJECT_NAME}"
  LANGUAGES CXX
  VERSION "${SKBUILD_PROJECT_VERSION}")

include(FetchContent)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(ABSL_PROPAGATE_CXX_STD ON)

# Remove support for scan deps, which is only useful when using C++ modules.
unset(CMAKE_CXX_SCANDEP_SOURCE)

# Do no build third party tests we are not using.
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
# Avoid creating a libcifpp cache directory with a redundant copy of the CCD.
set(CIFPP_CACHE_DIR "" CACHE PATH "" FORCE)


# Source archives instead of git tags so CMake doesn't do slow full git clones.
FetchContent_Declare(
  abseil-cpp
  URL https://github.com/abseil/abseil-cpp/archive/d7aaad83b488fd62bd51c81ecf16cd938532cc0a.tar.gz # 20240116.2
  URL_HASH SHA256=3a1f642ea21f45b220a01708d06b46042a3b750f31772b52371d6268f02d5bbc
  EXCLUDE_FROM_ALL)

FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/2e0815278cb899b20870a67ca8205996ef47e70f.tar.gz # v2.12.0
  URL_HASH SHA256=a50dfa85eda59b4c95180dd06baec89437ee1bc28ba92c5d9d1f853738ddd490
  EXCLUDE_FROM_ALL)

FetchContent_Declare(
  pybind11_abseil
  URL https://github.com/pybind/pybind11_abseil/archive/bddf30141f9fec8e577f515313caec45f559d319.tar.gz # HEAD @ 2024-08-07
  URL_HASH SHA256=a92e5faba0397eb4a6ba79317ab57de7891dc69ca8848cf2326abe28d5179c8f
  EXCLUDE_FROM_ALL)

# Transitive dependencies of libcifpp and dssp, declared also as source archives
# so CMake does not perform slow full git clones.
FetchContent_Declare(
  boost-rx
  URL https://github.com/boostorg/regex/archive/refs/tags/boost-1.83.0.tar.gz
  URL_HASH SHA256=1aa9cbe85e9a0077b20a3da3d1a48dbf2b9564d1c770e0cc5895354aff9f4616
  EXCLUDE_FROM_ALL)

FetchContent_Declare(
  my-eigen3
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
  URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
  EXCLUDE_FROM_ALL)

FetchContent_Declare(
  libmcfp
  URL https://github.com/mhekkel/libmcfp/archive/refs/tags/v1.3.1.tar.gz
  URL_HASH SHA256=28d298940d48610f45ea357c5b34c0d9e8c072cbb5e3582c27ce9a9038c1c400
  EXCLUDE_FROM_ALL)

FetchContent_Declare(
  cifpp
  URL https://github.com/pdb-redo/libcifpp/archive/ac98531a2fc8daf21131faa0c3d73766efa46180.tar.gz # v7.0.3
  URL_HASH SHA256=5fa208efa5aa7890b2a31ac1743ffb4dfc8bdc2605f590cc4639cc8922b18890
  # Don't `EXCLUDE_FROM_ALL` as necessary for build_data.
)

FetchContent_Declare(
  dssp
  URL https://github.com/PDB-REDO/dssp/archive/57560472b4260dc41f457706bc45fc6ef0bc0f10.tar.gz # v4.4.7
  URL_HASH SHA256=6cf3068527117e9a613728633a62f489af8890d5c30e7db919183b470e972258
  EXCLUDE_FROM_ALL)

FetchContent_MakeAvailable(pybind11 abseil-cpp pybind11_abseil cifpp dssp)

find_package(
  Python3
  COMPONENTS Interpreter Development NumPy
  REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(src/)

file(GLOB_RECURSE cpp_srcs src/alphafold3/*.cc)
list(FILTER cpp_srcs EXCLUDE REGEX ".*\(_test\|_main\|_benchmark\).cc$")

add_compile_definitions(NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

pybind11_add_module(cpp ${cpp_srcs})

target_link_libraries(
  cpp
  PRIVATE absl::check
          absl::flat_hash_map
          absl::node_hash_map
          absl::strings
          absl::status
          absl::statusor
          absl::log
          pybind11_abseil::absl_casters
          Python3::NumPy
          dssp::dssp
          cifpp::cifpp)

target_compile_definitions(cpp PRIVATE VERSION_INFO=${PROJECT_VERSION})
install(TARGETS cpp LIBRARY DESTINATION alphafold3)
install(
  FILES LICENSE
        OUTPUT_TERMS_OF_USE.md
        WEIGHTS_PROHIBITED_USE_POLICY.md
        WEIGHTS_TERMS_OF_USE.md
  DESTINATION alphafold3)
