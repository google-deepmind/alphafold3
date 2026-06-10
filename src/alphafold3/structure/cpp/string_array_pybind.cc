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

#include <Python.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

#include "numpy/arrayobject.h"
#include "numpy/ndarrayobject.h"
#include "numpy/ndarraytypes.h"
#include "numpy/npy_common.h"
#include "absl/algorithm/container.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "pybind11/cast.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"
#include "pybind11_abseil/absl_casters.h"

namespace {

namespace py = pybind11;

PyObject* RemapNumpyArrayObjects(PyObject* array, PyObject* mapping,
                                 bool inplace, PyObject* default_value) {
  import_array();
  if (!PyArray_Check(array)) {
    PyErr_SetString(PyExc_TypeError, "'array' must be a np.ndarray.");
    return nullptr;
  }
  if (!PyDict_Check(mapping)) {
    PyErr_SetString(PyExc_TypeError, "'mapping' must be a Python dict.");
    return nullptr;
  }

  PyArrayObject* array_obj = reinterpret_cast<PyArrayObject*>(array);
  if (PyArray_TYPE(array_obj) != NPY_OBJECT) {
    PyErr_SetString(PyExc_TypeError, "`array` must be an array of objects.");
    return nullptr;
  }

  if (inplace) {
    // We are returning original array so we need to increase the ref count.
    Py_INCREF(array);
  } else {
    // We are returning a fresh copy.
    array = PyArray_NewCopy(array_obj, NPY_CORDER);
    if (array == nullptr) {
      PyErr_SetString(PyExc_MemoryError, "Out of memory!");
      return nullptr;
    }
    array_obj = reinterpret_cast<PyArrayObject*>(array);
  }

  if (PyArray_SIZE(array_obj) == 0) {
    return array;
  }

  if (default_value == nullptr && PyDict_Size(mapping) == 0) {
    return array;
  }

  NpyIter* iter = NpyIter_New(
      array_obj, NPY_ITER_READWRITE | NPY_ITER_EXTERNAL_LOOP | NPY_ITER_REFS_OK,
      NPY_KEEPORDER, NPY_NO_CASTING, nullptr);
  if (iter == nullptr) {
    PyErr_SetString(PyExc_MemoryError, "Out of memory!");
    Py_XDECREF(array);
    return nullptr;
  }

  NpyIter_IterNextFunc* iter_next = NpyIter_GetIterNext(iter, nullptr);
  if (iter_next == nullptr) {
    NpyIter_Deallocate(iter);
    Py_XDECREF(array);
    PyErr_SetString(PyExc_MemoryError, "Out of memory!");
    return nullptr;
  }

  // Iterating arrays taken from:
  // https://numpy.org/doc/stable/reference/c-api/iterator.html
  char** data_pointer = NpyIter_GetDataPtrArray(iter);
  npy_intp* stride_pointer = NpyIter_GetInnerStrideArray(iter);
  npy_intp* inner_size_pointer = NpyIter_GetInnerLoopSizePtr(iter);
  do {
    char* data = *data_pointer;
    npy_intp stride = *stride_pointer;
    npy_intp count = *inner_size_pointer;
    for (size_t i = 0; i < count; ++i) {
      PyObject* entry;
      std::memcpy(&entry, data, sizeof(PyObject*));
      PyObject* result = PyDict_GetItem(mapping, entry);
      if (result != nullptr) {
        // Replace entry.
        Py_INCREF(result);
        Py_XDECREF(entry);
        std::memcpy(data, &result, sizeof(PyObject*));
      } else if (default_value != nullptr) {
        // Replace entry with a default value.
        Py_INCREF(default_value);
        Py_XDECREF(entry);
        std::memcpy(data, &default_value, sizeof(PyObject*));
      }
      data += stride;
    }
  } while (iter_next(iter));

  NpyIter_Deallocate(iter);
  return array;
}

// Convert 1D Numpy float array to a list of strings where each string has fixed
// number of decimal points. This is faster than Python list comprehension.
std::vector<std::string> FormatFloatArray(absl::Span<const float> values,
                                          int num_decimal_places) {
  std::vector<std::string> output;
  output.reserve(values.size());

  absl::c_transform(values, std::back_inserter(output),
                    [num_decimal_places](float value) {
                      return absl::StrFormat("%.*f", num_decimal_places, value);
                    });
  return output;
}

py::array_t<bool> IsIn(
    const py::array_t<PyObject*, py::array::c_style>& array,
    const absl::flat_hash_set<absl::string_view>& test_elements, bool invert) {
  const size_t num_elements = array.size();
  py::array_t<bool> output(num_elements);
  std::fill(output.mutable_data(), output.mutable_data() + output.size(),
            invert);

  // Shortcut: The output will be trivially always false if test_elements empty.
  if (test_elements.empty()) {
    return output;
  }

  for (size_t i = 0; i < num_elements; ++i) {
    // Compare the string values instead of comparing just object pointers.
    py::handle handle = array.data()[i];
    if (!PyUnicode_Check(handle.ptr()) && !PyBytes_Check(handle.ptr())) {
      continue;
    }
    if (test_elements.contains(py::cast<absl::string_view>(handle))) {
      output.mutable_data()[i] = !invert;
    }
  }
  if (array.ndim() > 1) {
    auto shape =
        std::vector<ptrdiff_t>(array.shape(), array.shape() + array.ndim());
    return output.reshape(shape);
  }
  return output;
}

py::array RemapMultipleArrays(
    const std::vector<py::array_t<PyObject*, py::array::c_style>>& arrays,
    const py::dict& mapping) {
  if (arrays.empty()) {
    return py::array_t<int64_t>(0);
  }
  size_t array_size = arrays[0].size();
  for (const auto& array : arrays) {
    if (array.size() != array_size) {
      throw py::value_error("All arrays must have the same length.");
    }
  }

  // Create a result buffer.
  auto result = py::array_t<int64_t>(array_size);
  absl::Span<int64_t> result_buffer(result.mutable_data(), array_size);
  std::vector<absl::Span<PyObject* const>> array_spans;
  array_spans.reserve(arrays.size());
  for (const auto& array : arrays) {
    array_spans.emplace_back(array.data(), array.size());
  }

  // Iterate over arrays and look up elements in the `py_dict`.
  bool fail = false;
  for (size_t i = 0; i < array_size; ++i) {
    PyObject* entry = PyTuple_New(arrays.size());
    if (entry == nullptr) {
      fail = true;
      break;
    }
    for (size_t j = 0; j < array_spans.size(); ++j) {
      PyTuple_SET_ITEM(entry, j, Py_NewRef(array_spans[j][i]));
    }
    PyObject* result = PyDict_GetItem(mapping.ptr(), entry);
    if (result != nullptr) {
      int64_t result_value = PyLong_AsLongLong(result);
      if (result_value == -1 && PyErr_Occurred()) {
        fail = true;
        Py_DECREF(entry);
        break;
      }
      if (result_value > std::numeric_limits<int64_t>::max() ||
          result_value < std::numeric_limits<int64_t>::lowest()) {
        PyErr_SetString(PyExc_OverflowError, "Result value too large.");
        fail = true;
        Py_DECREF(entry);
        break;
      }
      result_buffer[i] = result_value;
    } else {
      PyErr_Format(PyExc_KeyError, "%R", entry);
      fail = true;
      Py_DECREF(entry);
      break;
    }
    Py_DECREF(entry);
  }
  if (fail) {
    throw py::error_already_set();
  }
  return result;
}

constexpr char kRemapNumpyArrayObjects[] = R"(
Replace objects in NumPy array of objects using mapping.

Args:
  array: NumPy array with dtype=object.
  mapping: Dict mapping old values to new values.
  inplace: Bool (default False) whether to replace values inplace or to
    create a new array.
  default_value: If given, what value to map to if the mapping is missing
    for that particular item. If not given, such items are left unchanged.

Returns
  NumPy array of dtype object with values replaced according to mapping.
  If inplace is True the original array is modified inplace otherwise a
  new array is returned.
)";

constexpr char kFormatFloatArrayDoc[] = R"(
Converts float -> string array with given number of decimal places.
)";

constexpr char kIsInDoc[] = R"(
Computes whether each element is in test_elements.

Same use as np.isin, but much faster. If len(array) = n, len(test_elements) = m:
* This function has complexity O(n).
* np.isin with arrays of objects has complexity O(m*log(m) + n * log(m)).

Args:
  array: Input NumPy array with dtype=object.
  test_elements: The values against which to test each value of array.
  invert: If True, the values in the returned array are inverted, as if
    calculating `element not in test_elements`. Default is False.
    `isin(a, b, invert=True)` is equivalent to but faster than `~isin(a, b)`.

Returns
  A boolean array of the same shape as the input array. Each value `val` is:
  * `val in test_elements` if `invert=False`,
  * `val not in test_elements` if `invert=True`.
)";

constexpr char kRemapMultipleDoc[] = R"(
Maps keys from multiple aligned arrays to a single array.

Args:
  arrays: Numpy arrays of the same length. The tuple of aligned entries is used
    as key for the mapping.
  mapping: Dict mapping from tuples to integer values.

Returns
  NumPy array of dtype `int` with values looked up in mapping according to the
  tuple of aligned array entries as keys.
)";

}  // namespace

namespace alphafold3 {

void RegisterModuleStringArray(pybind11::module m) {
  m.def(
      "remap",
      [](py::object array, py::object mapping, bool inplace,
         py::object default_value) -> py::object {
        PyObject* result = RemapNumpyArrayObjects(array.ptr(), mapping.ptr(),
                                                  inplace, default_value.ptr());
        if (result == nullptr) {
          throw py::error_already_set();
        }
        return py::reinterpret_steal<py::object>(result);
      },
      py::return_value_policy::take_ownership, py::arg("array"),
      py::arg("mapping"), py::arg("inplace") = false, py::arg("default_value"),
      py::doc(kRemapNumpyArrayObjects + 1));
  m.def(
      "remap",
      [](py::object array, py::object mapping, bool inplace) -> py::object {
        PyObject* result = RemapNumpyArrayObjects(array.ptr(), mapping.ptr(),
                                                  inplace, nullptr);
        if (result == nullptr) {
          throw py::error_already_set();
        }
        return py::reinterpret_steal<py::object>(result);
      },
      py::return_value_policy::take_ownership, py::arg("array"),
      py::arg("mapping"), py::arg("inplace") = false,
      py::doc(kRemapNumpyArrayObjects + 1));
  m.def("format_float_array", &FormatFloatArray, py::arg("values"),
        py::arg("num_decimal_places"), py::doc(kFormatFloatArrayDoc + 1),
        py::call_guard<py::gil_scoped_release>());
  m.def("isin", &IsIn, py::arg("array"), py::arg("test_elements"),
        py::kw_only(), py::arg("invert") = false, py::doc(kIsInDoc + 1));
  m.def("remap_multiple", &RemapMultipleArrays, py::arg("arrays"),
        py::arg("mapping"), py::doc(kRemapMultipleDoc + 1));
}

}  // namespace alphafold3
