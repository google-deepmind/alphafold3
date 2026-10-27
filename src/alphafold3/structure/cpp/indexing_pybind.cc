// Copyright 2026 DeepMind Technologies Limited
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

#define NPY_TARGET_VERSION NPY_2_1_API_VERSION

#include "alphafold3/structure/cpp/indexing_pybind.h"

#include <Python.h>

#include <cstddef>
#include <cstdint>
#include <cstring>

#include "numpy/ndarrayobject.h"
#include "numpy/ndarraytypes.h"
#include "numpy/npy_common.h"
#include "absl/algorithm/container.h"
#include "absl/types/span.h"
#include "alphafold3/structure/cpp/key_lookup.h"
#include "pybind11/gil.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"

namespace {

namespace py = pybind11;

using ::alphafold3::ResolveKeys;

// Returns a native-endian, aligned, C-contiguous int64 view or copy of `array`.
py::array_t<int64_t> AsAlignedCContiguousInt64(const py::array& array) {
  auto* obj = reinterpret_cast<PyArrayObject*>(array.ptr());
  if (PyArray_ISCARRAY_RO(obj) && PyArray_ISSIGNED(obj) &&
      PyArray_ITEMSIZE(obj) == sizeof(int64_t)) {
    return py::reinterpret_borrow<py::array_t<int64_t>>(array);
  }
  PyObject* copy =
      PyArray_FROM_OTF(array.ptr(), NPY_INT64, NPY_ARRAY_CARRAY_RO);
  if (copy == nullptr) {
    throw py::error_already_set();
  }
  return py::reinterpret_steal<py::array_t<int64_t>>(copy);
}

absl::Span<const int64_t> AsSpan(const py::array_t<int64_t>& a) {
  return {a.data(), static_cast<size_t>(a.size())};
}

[[noreturn]] void ThrowMissingKeys() {
  throw py::index_error("Keys missing from table.");
}

void ValidateKeyArrays(const py::array& keys, const py::array& table_keys) {
  if (!PyArray_ISSIGNED(reinterpret_cast<PyArrayObject*>(keys.ptr())) ||
      !PyArray_ISSIGNED(reinterpret_cast<PyArrayObject*>(table_keys.ptr())) ||
      keys.ndim() != 1 || table_keys.ndim() != 1) {
    throw py::index_error("Keys must be 1D signed-integer arrays.");
  }
}

bool IsIdentityKeys(absl::Span<const int64_t> table_span) {
  if (table_span.front() != 0 ||
      table_span.back() != static_cast<int64_t>(table_span.size() - 1)) {
    return false;
  }
  for (size_t i = 1; i + 1 < table_span.size(); ++i) {
    if (table_span[i] != static_cast<int64_t>(i)) {
      return false;
    }
  }
  return true;
}

py::array_t<int64_t> IndexOf(const py::array& keys,
                             const py::array& table_keys) {
  ValidateKeyArrays(keys, table_keys);
  if (keys.size() > 0 && table_keys.size() == 0) {
    ThrowMissingKeys();
  }
  if (keys.size() == 0) {
    return py::array_t<int64_t>(0);
  }
  auto contig_table = AsAlignedCContiguousInt64(table_keys);
  auto contig_keys = AsAlignedCContiguousInt64(keys);
  const auto table_span = AsSpan(contig_table);
  const auto key_span = AsSpan(contig_keys);
  const bool can_borrow = contig_keys.ptr() != keys.ptr() || !keys.writeable();
  bool all_present = false;
  if (can_borrow) {
    bool is_identity = false;
    {
      py::gil_scoped_release release;
      is_identity = IsIdentityKeys(table_span);
      if (is_identity) {
        const uint64_t num_rows = table_span.size();
        all_present = absl::c_all_of(key_span, [num_rows](int64_t k) {
          return static_cast<uint64_t>(k) < num_rows;
        });
      }
    }
    if (is_identity) {
      if (!all_present) {
        ThrowMissingKeys();
      }
      return contig_keys;
    }
  }
  py::array_t<int64_t> indices(keys.size());
  auto indices_span =
      absl::MakeSpan(indices.mutable_data(), static_cast<size_t>(keys.size()));
  {
    py::gil_scoped_release release;
    all_present = ResolveKeys(key_span, table_span, indices_span);
  }
  if (!all_present) {
    ThrowMissingKeys();
  }
  return indices;
}

py::object NonzeroIndices(const py::array& mask) {
  static_assert(sizeof(bool) == 1);
  auto* obj = reinterpret_cast<PyArrayObject*>(mask.ptr());
  if (!PyArray_ISBOOL(obj) || mask.ndim() != 1) {
    throw py::index_error("Mask must be a 1D bool array.");
  }
  auto contig = py::array_t<bool, py::array::c_style>::ensure(mask);
  static_assert(sizeof(bool) == 1);
  const auto span =
      absl::MakeConstSpan(reinterpret_cast<const uint8_t*>(contig.data()),
                          static_cast<size_t>(contig.size()));
  size_t nnz = 0;
  {
    py::gil_scoped_release release;
    nnz = static_cast<size_t>(
        absl::c_count_if(span, [](uint8_t b) { return b != 0; }));
  }
  if (nnz == span.size()) {
    return py::none();
  }
  py::array_t<int64_t> out(static_cast<py::ssize_t>(nnz));
  int64_t* dst = out.mutable_data();
  int64_t* const dst_end = dst + nnz;
  {
    py::gil_scoped_release release;
    for (size_t i = 0, n = span.size(); i < n;) {
      if (i + sizeof(uint64_t) <= n) {
        uint64_t word = 0;
        std::memcpy(&word, span.data() + i, sizeof(word));
        if (word == 0) {
          i += sizeof(uint64_t);
          continue;
        }
      }
      if (span[i] != 0) {
        *dst++ = static_cast<int64_t>(i);
      }
      ++i;
    }
  }
  if (dst != dst_end) {
    throw py::value_error("Mask modified concurrently during nonzero_indices.");
  }
  return out;
}

}  // namespace

namespace alphafold3 {

constexpr char kIndexOfDoc[] =
    R"(Maps 1D signed-integer keys to their 0-based positions in 1D table_keys.

Args:
  keys: 1D signed-integer NumPy array of query keys.
  table_keys: 1D signed-integer NumPy array of table keys (if duplicate keys
    exist, the last matching index is returned).

Returns:
  A 1D int64 NumPy array with the same length as `keys` containing the 0-based
  index in `table_keys` for each element of `keys` (may alias `keys` when it
  is a read-only int64 array and `table_keys` is `0..m - 1`).

Raises:
  IndexError: If `keys` or `table_keys` is not a 1D signed-integer array, or
    any element of `keys` is absent from `table_keys`.
)";

constexpr char kNonzeroIndicesDoc[] =
    R"(Returns int64 indices of True elements in a 1D bool mask, or None if all True.)";

void RegisterModuleIndexing(pybind11::module m) {
  if (_import_array() < 0) {
    throw py::error_already_set();
  }
  m.def("index_of", &IndexOf, py::arg("keys").noconvert(),
        py::arg("table_keys").noconvert(), py::doc(kIndexOfDoc));
  m.def("nonzero_indices", &NonzeroIndices, py::arg("mask").noconvert(),
        py::doc(kNonzeroIndicesDoc));
}

}  // namespace alphafold3
