# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
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

FetchContent_Declare(
  abseil-cpp
  GIT_REPOSITORY https://github.com/abseil/abseil-cpp
  GIT_TAG d7aaad83b488fd62bd51c81ecf16cd938532cc0a # 20240116.2
  EXCLUDE_FROM_ALL)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11
  GIT_TAG 2e0815278cb899b20870a67ca8205996ef47e70f # v2.12.0
  EXCLUDE_FROM_ALL)

FetchContent_Declare(
  pybind11_abseil
  GIT_REPOSITORY https://github.com/pybind/pybind11_abseil
  GIT_TAG bddf30141f9fec8e577f515313caec45f559d319 # HEAD @ 2024-08-07
  EXCLUDE_FROM_ALL)

FetchContent_Declare(
  cifpp
  GIT_REPOSITORY https://github.com/pdb-redo/libcifpp
  GIT_TAG ac98531a2fc8daf21131faa0c3d73766efa46180 # v7.0.3
  # Don't `EXCLUDE_FROM_ALL` as necessary for build_data.
)

FetchContent_Declare(
  dssp
  GIT_REPOSITORY https://github.com/PDB-REDO/dssp
  GIT_TAG 57560472b4260dc41f457706bc45fc6ef0bc0f10 # v4.4.7
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
