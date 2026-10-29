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

#define NPY_TARGET_VERSION NPY_2_1_API_VERSION

#include "alphafold3/structure/cpp/membership_pybind.h"

#include <Python.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>

#include "numpy/ndarrayobject.h"
#include "numpy/ndarraytypes.h"
#include "numpy/npy_common.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "alphafold3/structure/cpp/key_lookup.h"
#include "pybind11/cast.h"
#include "pybind11/gil.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"
#include "pybind11_abseil/absl_casters.h"

namespace {

namespace py = pybind11;

using ::alphafold3::ContainsKeys;

// Returns a native-endian, aligned, C-contiguous view or copy of `array`.
template <typename T>
py::array_t<T> AsAlignedCContiguous(const py::array& array) {
  constexpr int kTargetTypeNum =
      std::is_same_v<T, PyObject*>
          ? NPY_OBJECT
          : (sizeof(T) == sizeof(int64_t) ? NPY_INT64 : NPY_INT32);
  auto* arr_obj = reinterpret_cast<PyArrayObject*>(array.ptr());
  const bool matches_type = std::is_same_v<T, PyObject*>
                                ? PyArray_ISOBJECT(arr_obj)
                                : (PyArray_ISSIGNED(arr_obj) &&
                                   PyArray_ITEMSIZE(arr_obj) == sizeof(T));
  if (PyArray_ISCARRAY_RO(arr_obj) && matches_type) {
    return py::reinterpret_borrow<py::array_t<T>>(array);
  }
  PyObject* obj =
      PyArray_FROM_OTF(array.ptr(), kTargetTypeNum, NPY_ARRAY_CARRAY_RO);
  if (obj == nullptr) {
    throw py::error_already_set();
  }
  return py::reinterpret_steal<py::array_t<T>>(obj);
}

py::array_t<bool> NewBoolCArray(const py::array& shape_src) {
  PyObject* obj =
      PyArray_SimpleNew(shape_src.ndim(), shape_src.shape(), NPY_BOOL);
  if (obj == nullptr) {
    throw py::error_already_set();
  }
  return py::reinterpret_steal<py::array_t<bool>>(obj);
}

py::array_t<bool> MakeConstantBoolArray(const py::array& shape_src,
                                        bool value) {
  auto result = NewBoolCArray(shape_src);
  std::fill_n(result.mutable_data(), result.size(), value);
  return result;
}

py::array_t<bool> IsInObjectSet(const py::array& array, PyObject* set_ptr,
                                bool invert) {
  auto contig = AsAlignedCContiguous<PyObject*>(array);
  auto result = NewBoolCArray(contig);
  bool* out = result.mutable_data();
  PyObject* const* vals = contig.data();
  py::object prev_val;
  bool prev_out = invert;
  for (size_t i = 0, n = contig.size(); i < n; ++i) {
    if (PyObject* val = vals[i]; val != prev_val.ptr()) {
      prev_val = py::reinterpret_borrow<py::object>(val);
      if (val == nullptr) {
        prev_out = invert;
      } else {
        const int contains = PySet_Contains(set_ptr, val);
        if (contains < 0) {
          throw py::error_already_set();
        }
        prev_out = invert != (contains != 0);
      }
    }
    out[i] = prev_out;
  }
  return result;
}

// Skips out-of-range integers while keeping `py::cast<IntT>` non-int errors.
template <typename IntT>
absl::flat_hash_set<IntT> CastIntSet(const py::iterable& test_elements) {
  absl::flat_hash_set<IntT> result;
  result.reserve(py::len_hint(test_elements));
  for (py::handle elem : test_elements) {
    if (!PyIndex_Check(elem.ptr())) {
      py::cast<IntT>(elem);
    }
    int overflow = 0;
    const int64_t val = PyLong_AsLongLongAndOverflow(elem.ptr(), &overflow);
    if (PyErr_Occurred()) {
      throw py::error_already_set();
    }
    if (overflow == 0 && val >= std::numeric_limits<IntT>::min() &&
        val <= std::numeric_limits<IntT>::max()) {
      result.insert(static_cast<IntT>(val));
    }
  }
  return result;
}

template <typename IntT>
py::array_t<bool> IsInIntSet(const py::array& array,
                             const py::anyset& test_elements, bool invert) {
  auto contig = AsAlignedCContiguous<IntT>(array);
  absl::flat_hash_set<IntT> test_set =
      CastIntSet<IntT>(py::reinterpret_borrow<py::iterable>(test_elements));
  auto result = NewBoolCArray(contig);
  bool* out = result.mutable_data();
  const IntT* vals = contig.data();
  const size_t n = contig.size();
  {
    py::gil_scoped_release release;
    for (size_t i = 0; i < n; ++i) {
      out[i] = invert != test_set.contains(vals[i]);
    }
  }
  return result;
}

template <typename IntT>
py::array_t<bool> IsInTypedIntArray(const py::array& array,
                                    const py::array& test_elements,
                                    bool invert) {
  auto contig = AsAlignedCContiguous<IntT>(array);
  auto contig_test = AsAlignedCContiguous<IntT>(test_elements);
  auto result = NewBoolCArray(array);
  {
    py::gil_scoped_release release;
    ContainsKeys<IntT>(absl::MakeConstSpan(contig),
                       absl::MakeConstSpan(contig_test), invert,
                       absl::MakeSpan(result.mutable_data(), result.size()));
  }
  return result;
}

py::array_t<bool> IsInIntArray(const py::array& array,
                               const py::array& test_elements, bool invert) {
  auto* arr_obj = reinterpret_cast<PyArrayObject*>(array.ptr());
  auto* test_obj = reinterpret_cast<PyArrayObject*>(test_elements.ptr());
  if (!PyArray_ISSIGNED(arr_obj) || !PyArray_ISSIGNED(test_obj)) {
    throw py::type_error("Incompatible array and test_elements dtypes.");
  }
  if (array.size() == 0 || test_elements.size() == 0) {
    return MakeConstantBoolArray(array, invert);
  }
  if (PyArray_ITEMSIZE(arr_obj) <= sizeof(int32_t) &&
      PyArray_ITEMSIZE(test_obj) <= sizeof(int32_t)) {
    return IsInTypedIntArray<int32_t>(array, test_elements, invert);
  }
  return IsInTypedIntArray<int64_t>(array, test_elements, invert);
}

absl::string_view ToStringView(const char* buf, size_t size) {
  return size == 0 ? absl::string_view() : absl::string_view(buf, size);
}

// `iter` is a two-operand iterator created by `NpyIter_MultiNew`, where
// operand 0 is the read-only NPY_VSTRING input array and operand 1 is the
// NPY_BOOL output array allocated by the iterator. `descr` must be the
// descriptor of operand 0, and owns the string allocator holding the packed
// string data.
void PopulateStringIsInResult(
    NpyIter* iter, PyArray_Descr* descr,
    const absl::flat_hash_set<std::string>& test_elements_set, bool invert) {
  if (NpyIter_GetIterSize(iter) == 0) {
    return;
  }
  if (test_elements_set.empty()) {
    if (PyArray_FillWithScalar(NpyIter_GetOperandArray(iter)[1],
                               py::bool_(invert).ptr()) < 0) {
      throw py::error_already_set();
    }
    return;
  }

  // Stored in a local variable to avoid reinitializing it every time.
  NpyIter_IterNextFunc* iter_next = NpyIter_GetIterNext(iter, nullptr);
  if (iter_next == nullptr) {
    // NpyIter_GetIterNext sets the python error.
    throw py::error_already_set();
  }

  // Iterating arrays taken from:
  // https://numpy.org/doc/stable/reference/c-api/iterator.html
  char** data_pointer = NpyIter_GetDataPtrArray(iter);
  npy_intp* stride_pointer = NpyIter_GetInnerStrideArray(iter);
  npy_intp* inner_size_pointer = NpyIter_GetInnerLoopSizePtr(iter);

  int load_status = 0;
  {
    py::gil_scoped_release release;
    PyArray_StringDTypeObject* descr_as_array =
        reinterpret_cast<PyArray_StringDTypeObject*>(descr);
    const bool null_in_set = descr_as_array->has_string_na &&
                             test_elements_set.contains(ToStringView(
                                 descr_as_array->default_string.buf,
                                 descr_as_array->default_string.size));
    std::unique_ptr<npy_string_allocator,
                    decltype(&NpyString_release_allocator)>
        allocator(NpyString_acquire_allocator(descr_as_array),
                  &NpyString_release_allocator);
    do {
      const char* string_data = data_pointer[0];
      char* result_data = data_pointer[1];
      const npy_intp string_stride = stride_pointer[0];
      const npy_intp result_stride = stride_pointer[1];
      const npy_intp count = *inner_size_pointer;
      for (npy_intp i = 0; i < count; ++i) {
        const auto* current_packed_string =
            reinterpret_cast<const npy_packed_static_string*>(string_data);
        npy_static_string str = {0, nullptr};
        load_status =
            NpyString_load(allocator.get(), current_packed_string, &str);
        if (load_status < 0) {
          break;
        }
        const bool in_set =
            load_status == 0
                ? test_elements_set.contains(ToStringView(str.buf, str.size))
                : (load_status == 1 && null_in_set);
        *reinterpret_cast<npy_bool*>(result_data) = invert != in_set;
        string_data += string_stride;
        result_data += result_stride;
      }
    } while (load_status >= 0 && iter_next(iter));
  }
  if (load_status < 0) {
    throw py::value_error("Failed to load StringDType string.");
  }
}

py::array_t<bool> IsIn(py::handle array_arg, const py::anyset& test_elements,
                       bool invert) {
  py::array array = PyArray_CheckExact(array_arg.ptr())
                        ? py::reinterpret_borrow<py::array>(array_arg)
                        : py::array::ensure(array_arg);
  if (!array) {
    throw py::type_error("Incompatible array argument.");
  }

  auto* array_obj = reinterpret_cast<PyArrayObject*>(array.ptr());
  PyArray_Descr* descr = PyArray_DESCR(array_obj);

  if (descr->type_num == NPY_VSTRING) {
    absl::flat_hash_set<std::string> test_elements_set;
    test_elements_set.reserve(PySet_GET_SIZE(test_elements.ptr()));
    for (py::handle elem : test_elements) {
      if (PyUnicode_Check(elem.ptr())) {
        Py_ssize_t size = 0;
        const char* data = PyUnicode_AsUTF8AndSize(elem.ptr(), &size);
        if (data == nullptr) {
          throw py::error_already_set();
        }
        test_elements_set.emplace(data, static_cast<size_t>(size));
      }
    }

    PyArrayObject* op[2] = {array_obj, nullptr};
    npy_uint32 op_flags[2] = {
        NPY_ITER_READONLY,
        NPY_ITER_WRITEONLY | NPY_ITER_ALLOCATE,
    };
    constexpr npy_uint32 kFlags =
        NPY_ITER_ZEROSIZE_OK | NPY_ITER_EXTERNAL_LOOP | NPY_ITER_REFS_OK;

    py::dtype bool_dtype = py::dtype::of<bool>();
    PyArray_Descr* op_dtypes[2] = {
        nullptr, reinterpret_cast<PyArray_Descr*>(bool_dtype.ptr())};

    std::unique_ptr<NpyIter, decltype(&NpyIter_Deallocate)> iter(
        NpyIter_MultiNew(
            /*nop=*/2, op, kFlags,
            /*order=*/NPY_KEEPORDER,
            /*casting=*/NPY_NO_CASTING, op_flags, op_dtypes),
        &NpyIter_Deallocate);

    if (iter == nullptr) {
      // NpyIter_MultiNew sets the python error.
      throw py::error_already_set();
    }

    PopulateStringIsInResult(iter.get(), descr, test_elements_set, invert);

    py::array_t<bool> result = py::reinterpret_borrow<py::array_t<bool>>(
        reinterpret_cast<PyObject*>(NpyIter_GetOperandArray(iter.get())[1]));

    // Circumvent iterator deallocation in the successful case so that we can
    // correctly check deallocation status.
    if (NpyIter_Deallocate(iter.release()) != NPY_SUCCEED) {
      throw py::error_already_set();
    }

    return result;
  }

  if (test_elements.empty()) {
    return MakeConstantBoolArray(array, invert);
  }
  const npy_intp elsize = PyArray_ITEMSIZE(array_obj);
  // Reuse `int32_t`/`int64_t` kernels (`uint64` falls back to `IsInObjectSet`).
  if (PyDataType_ISSIGNED(descr) && elsize <= sizeof(int32_t)) {
    return IsInIntSet<int32_t>(array, test_elements, invert);
  }
  if ((PyDataType_ISSIGNED(descr) && elsize == sizeof(int64_t)) ||
      (PyDataType_ISUNSIGNED(descr) && elsize <= sizeof(uint32_t))) {
    return IsInIntSet<int64_t>(array, test_elements, invert);
  }
  return IsInObjectSet(array, test_elements.ptr(), invert);
}

constexpr char kIsInDoc[] =
    R"(Computes whether each element is in test_elements.

Same use as np.isin, but much faster. If len(array) = n, len(test_elements) = m:
* This function has complexity O(n).
* np.isin with kind='sort' has complexity O(m * log(m) + n * log(m)).

Args:
  array: Input NumPy array. The type of the member elements must be compatible
    with the type of the elements in test_elements.
  test_elements: A Python set/frozenset, or a signed-integer NumPy array when
    `array` is a signed-integer array.
  invert: If True, the values in the returned array are inverted, as if
    calculating `element not in test_elements`. Default is False. `isin(a, b,
    invert=True)` is equivalent to but faster than `~isin(a, b)`.

Returns
  A boolean array of the same shape as the input array. Each value `val` is:
  * `val in test_elements` if `invert=False`,
  * `val not in test_elements` if `invert=True`.
)";

}  // namespace

namespace alphafold3 {

void RegisterModuleMembership(pybind11::module m) {
  if (_import_array() < 0) {
    throw py::error_already_set();
  }
  m.def("isin", &IsInIntArray, py::arg("array").noconvert(),
        py::arg("test_elements").noconvert(), py::kw_only(),
        py::arg("invert") = false, py::doc(kIsInDoc));
  m.def("isin", &IsIn, py::arg("array"), py::arg("test_elements"),
        py::kw_only(), py::arg("invert") = false, py::doc(kIsInDoc));
}

}  // namespace alphafold3
