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

#include <Python.h>

#include "alphafold3/structure/cpp/membership_pybind.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>

#include "numpy/ndarrayobject.h"
#include "numpy/ndarraytypes.h"
#include "numpy/npy_common.h"
#include "absl/algorithm/container.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "pybind11/cast.h"
#include "pybind11/gil.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"
#include "pybind11_abseil/absl_casters.h"

namespace {

namespace py = pybind11;

void ApplyIsIn(absl::Span<bool> result, absl::Span<PyObject* const> values,
               py::anyset test_elements, bool invert) {
  for (size_t i = 0, e = values.size(); i < e; ++i) {
    if (values[i] == nullptr) {
      result[i] = invert;
      continue;
    }
    py::handle test_element = py::cast<py::handle>(values[i]);
    result[i] = invert != test_elements.contains(test_element);
  }
}

// Skips out-of-range integers while keeping `py::cast<Number>` non-int errors.
template <typename Number>
absl::flat_hash_set<Number> CastIntSet(const py::iterable& test_elements) {
  absl::flat_hash_set<Number> result;
  result.reserve(py::len_hint(test_elements));
  for (py::handle elem : test_elements) {
    if (!PyIndex_Check(elem.ptr())) {
      py::cast<Number>(elem);
    }
    int overflow = 0;
    const int64_t val = PyLong_AsLongLongAndOverflow(elem.ptr(), &overflow);
    if (PyErr_Occurred()) {
      throw py::error_already_set();
    }
    if (overflow == 0 && val >= std::numeric_limits<Number>::min() &&
        val <= std::numeric_limits<Number>::max()) {
      result.insert(static_cast<Number>(val));
    }
  }
  return result;
}

template <typename Number>
void ApplyIsIn(absl::Span<bool> result, absl::Span<const Number> values,
               py::anyset test_elements, bool invert) {
  auto test_elements_set =
      CastIntSet<Number>(py::reinterpret_borrow<py::iterable>(test_elements));
  py::gil_scoped_release release;
  for (size_t i = 0, e = values.size(); i < e; ++i) {
    result[i] = invert != test_elements_set.contains(values[i]);
  }
}

template <typename T>
py::array_t<bool> IsIn(const py::array_t<T, py::array::c_style>& array,
                       py::anyset test_elements, bool invert) {
  py::array_t<bool, py::array::c_style> result(
      absl::MakeConstSpan(array.shape(), array.ndim()));
  const size_t size = array.size();
  auto result_span = absl::MakeSpan(result.mutable_data(), size);
  if (test_elements.attr("__len__")().cast<ssize_t>() == 0) {
    absl::c_fill(result_span, invert);
  } else {
    auto array_span = absl::MakeConstSpan(array.data(), size);
    ApplyIsIn(result_span, array_span, std::move(test_elements), invert);
  }
  return result;
}

// `iter` is a two-operand iterator created by `NpyIter_MultiNew`, where
// operand 0 is the read-only NPY_VSTRING input array and operand 1 is the
// NPY_BOOL output array allocated by the iterator. `descr` must be the
// descriptor of operand 0, and owns the string allocator holding the packed
// string data.
void PopulateStringIsInResult(
    NpyIter* iter, PyArray_Descr* descr,
    const absl::flat_hash_set<absl::string_view>& test_elements_set,
    bool invert) {
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

  PyArray_StringDTypeObject* descr_as_array =
      reinterpret_cast<PyArray_StringDTypeObject*>(descr);
  std::unique_ptr<npy_string_allocator, decltype(&NpyString_release_allocator)>
      allocator(NpyString_acquire_allocator(descr_as_array),
                &NpyString_release_allocator);
  if (allocator == nullptr) {
    throw py::error_already_set();
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
  do {
    const char* string_data = data_pointer[0];
    char* result_data = data_pointer[1];
    npy_intp string_stride = stride_pointer[0];
    npy_intp result_stride = stride_pointer[1];
    npy_intp count = *inner_size_pointer;
    for (npy_intp i = 0; i < count; ++i) {
      // Access StringDType array element.
      const auto* current_packed_string =
          reinterpret_cast<const npy_packed_static_string*>(string_data);
      npy_static_string str = {0, nullptr};
      switch (NpyString_load(allocator.get(), current_packed_string, &str)) {
        case 0:  // The string is not null.
          *reinterpret_cast<npy_bool*>(result_data) =
              invert !=
              test_elements_set.contains(absl::string_view(str.buf, str.size));
          break;
        case 1:  // String is null.
          *reinterpret_cast<npy_bool*>(result_data) = invert;
          break;
        default:
          throw py::error_already_set();
      }
      string_data += string_stride;
      result_data += result_stride;
    }
  } while (iter_next(iter));
}

py::array_t<bool> isin(const py::array& array, py::anyset test_elements,
                       bool invert) {
  // Ensure the array is of the correct StringDType
  PyArrayObject* array_obj = reinterpret_cast<PyArrayObject*>(array.ptr());
  PyArray_Descr* descr = PyArray_DESCR(array_obj);

  if (descr->type_num != NPY_VSTRING) {
    const npy_intp elsize = PyDataType_ELSIZE(descr);
    // Reuse `int32_t`/`int64_t` kernels (`uint64` falls back to `PyObject*`).
    if (PyDataType_ISSIGNED(descr) && elsize <= sizeof(int32_t)) {
      return IsIn(py::cast<py::array_t<int32_t, py::array::c_style>>(array),
                  std::move(test_elements), invert);
    }
    if ((PyDataType_ISSIGNED(descr) && elsize == sizeof(int64_t)) ||
        (PyDataType_ISUNSIGNED(descr) && elsize <= sizeof(uint32_t))) {
      return IsIn(py::cast<py::array_t<int64_t, py::array::c_style>>(array),
                  std::move(test_elements), invert);
    }
    // Fall back to the object-based IsIn, for backwards compatibility.
    return IsIn(py::cast<py::array_t<PyObject*, py::array::c_style>>(array),
                std::move(test_elements), invert);
  }

  absl::flat_hash_set<absl::string_view> test_elements_set =
      py::cast<absl::flat_hash_set<absl::string_view>>(test_elements);

  PyArrayObject* op[2] = {array_obj, nullptr};
  npy_uint32 op_flags[2] = {
      NPY_ITER_READONLY,
      NPY_ITER_WRITEONLY | NPY_ITER_ALLOCATE,
  };
  npy_uint32 flags =
      NPY_ITER_ZEROSIZE_OK | NPY_ITER_EXTERNAL_LOOP | NPY_ITER_REFS_OK;

  py::dtype bool_dtype = py::dtype::of<bool>();
  PyArray_Descr* op_dtypes[2] = {
      nullptr, reinterpret_cast<PyArray_Descr*>(bool_dtype.ptr())};

  std::unique_ptr<NpyIter, decltype(&NpyIter_Deallocate)> iter(
      NpyIter_MultiNew(
          /*nop=*/2, op, flags,
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

constexpr char kIsInDoc[] = R"(
Computes whether each element is in test_elements.

Same use as np.isin, but much faster. If len(array) = n, len(test_elements) = m:
* This function has complexity O(n).
* np.isin with kind='sort' has complexity O(m * log(m) + n * log(m)).

Args:
  array: Input NumPy array. The type of the member elements must be compatible
    with the type of the elements in test_elements.
  test_elements: The values against which to test each value of array.
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
  m.def("isin", &IsIn<int64_t>, py::arg("array"), py::arg("test_elements"),
        py::kw_only(), py::arg("invert") = false, py::doc(kIsInDoc + 1));
  m.def("isin", &IsIn<int32_t>, py::arg("array"), py::arg("test_elements"),
        py::kw_only(), py::arg("invert") = false, py::doc(kIsInDoc + 1));
  m.def("isin", &IsIn<PyObject*>, py::arg("array"), py::arg("test_elements"),
        py::kw_only(), py::arg("invert") = false, py::doc(kIsInDoc + 1));
  m.def("isin", &isin, py::arg("array"), py::arg("test_elements"),
        py::kw_only(), py::arg("invert") = false, py::doc(kIsInDoc + 1));
}

}  // namespace alphafold3
