// Copyright 2024 DeepMind Technologies Limited
//
// AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
// this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// To request access to the AlphaFold 3 model parameters, follow the process set
// out at https://github.com/google-deepmind/alphafold3. You may only use these
// if received directly from Google. Use is subject to terms of use available at
// https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

#ifndef ALPHAFOLD3_SRC_ALPHAFOLD3_MODEL_JSON_SERIALIZE_PYBIND_H_
#define ALPHAFOLD3_SRC_ALPHAFOLD3_MODEL_JSON_SERIALIZE_PYBIND_H_

#include "pybind11/pybind11.h"

namespace alphafold3 {

void RegisterModuleJsonSerialize(pybind11::module m);

}

#endif  // ALPHAFOLD3_SRC_ALPHAFOLD3_MODEL_JSON_SERIALIZE_PYBIND_H_
