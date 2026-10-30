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
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "numpy/ndarrayobject.h"
#include "numpy/ndarraytypes.h"
#include "numpy/npy_common.h"
#include "absl/base/no_destructor.h"
#include "absl/container/flat_hash_map.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "alphafold3/parsers/cpp/cif_dict_lib.h"
#include "pybind11/attr.h"
#include "pybind11/cast.h"
#include "pybind11/gil.h"
#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"
#include "pybind11/stl.h"

namespace alphafold3 {
namespace {
namespace py = pybind11;

template <typename Item, typename ForEach>
bool GatherArray(size_t num_dims, const npy_intp* shape_array,
                 const npy_intp* stride_array, const char* data,
                 absl::Span<const std::string> values,
                 const ForEach& for_each_cb, std::string& index_error) {
  if (num_dims <= 1) {
    const npy_intp shape = num_dims == 1 ? shape_array[0] : 1;
    const npy_intp stride = num_dims == 1 ? stride_array[0] : 0;
    for (npy_intp i = 0; i < shape; ++i) {
      Item index;
      std::memcpy(&index, data + stride * i, sizeof(Item));
      if (index < 0 || static_cast<size_t>(index) >= values.size()) {
        index_error = absl::StrCat("index ", index,
                                   " is out of bounds for column with size ",
                                   values.size());
        return false;
      }
      if (!for_each_cb(values[index])) {
        return false;
      }
    }
    return true;
  }
  const npy_intp shape = shape_array[0];
  const npy_intp stride = stride_array[0];
  for (npy_intp i = 0; i < shape; ++i) {
    if (!GatherArray<Item>(num_dims - 1, shape_array + 1, stride_array + 1,
                           data + stride * i, values, for_each_cb,
                           index_error)) {
      return false;
    }
  }
  return true;
}

template <bool kReleaseGil = false, typename Size, typename ForEach>
bool Gather(PyObject* gather, absl::Span<const std::string> values,
            Size&& size_cb, const ForEach& for_each_cb) {
  auto run_loop = [](auto&& fn) {
    if constexpr (kReleaseGil) {
      py::gil_scoped_release gil_release;
      return fn();
    } else {
      return fn();
    }
  };
  if (gather == Py_None) {
    npy_intp dim = static_cast<npy_intp>(values.size());
    if (!size_cb(absl::MakeSpan(&dim, 1))) {
      return false;
    }
    return run_loop([&] {
      for (const std::string& v : values) {
        if (!for_each_cb(v)) {
          return false;
        }
      }
      return true;
    });
  }
  if (PySlice_Check(gather)) {
    Py_ssize_t start, stop, step, slice_length;
    if (PySlice_GetIndicesEx(gather, values.size(), &start, &stop, &step,
                             &slice_length) != 0) {
      return false;
    }
    npy_intp dim = static_cast<npy_intp>(slice_length);
    if (!size_cb(absl::MakeSpan(&dim, 1))) {
      return false;
    }
    return run_loop([&] {
      for (Py_ssize_t i = 0; i < slice_length; ++i) {
        if (!for_each_cb(values[start + i * step])) {
          return false;
        }
      }
      return true;
    });
  }
  if (PyArray_Check(gather)) {
    PyArrayObject* gather_array = reinterpret_cast<PyArrayObject*>(gather);
    const int ndim = PyArray_NDIM(gather_array);
    const npy_intp* shape = PyArray_DIMS(gather_array);
    const npy_intp* strides = PyArray_STRIDES(gather_array);
    const char* bytes = PyArray_BYTES(gather_array);
    auto gather_typed = [&](auto item_tag) {
      using Item = decltype(item_tag);
      if (!size_cb(absl::MakeSpan(shape, ndim))) {
        return false;
      }
      std::string index_error;
      const bool ok = run_loop([&] {
        return GatherArray<Item>(ndim, shape, strides, bytes, values,
                                 for_each_cb, index_error);
      });
      if (!ok && !index_error.empty()) {
        PyErr_SetString(PyExc_IndexError, index_error.c_str());
      }
      return ok;
    };
    switch (PyArray_TYPE(gather_array)) {
      case NPY_INT16:
        return gather_typed(int16_t{});
      case NPY_UINT16:
        return gather_typed(uint16_t{});
      case NPY_INT32:
        return gather_typed(int32_t{});
      case NPY_UINT32:
        return gather_typed(uint32_t{});
      case NPY_INT64:
        return gather_typed(int64_t{});
      case NPY_UINT64:
        return gather_typed(uint64_t{});
      default:
        PyErr_SetString(PyExc_TypeError, "Unsupported NumPy array type.");
        return false;
    }
  }

  PyErr_Format(PyExc_TypeError, "Invalid gather %R", gather);
  return false;
}

// Creates a NumPy array of objects of given strings. Reusing duplicates where
// possible.
PyObject* ConvertStrings(PyObject* gather, PyArray_Descr* type,
                         absl::Span<const std::string> values) {
  py::object type_owner =
      py::reinterpret_steal<py::object>(reinterpret_cast<PyObject*>(type));
  absl::flat_hash_map<absl::string_view, PyObject*> existing;
  existing.reserve(std::min<size_t>(values.size(), 64 * 1024));
  absl::string_view last_value;
  PyObject* last_obj = nullptr;

  PyObject* ret = nullptr;
  PyObject** dst = nullptr;
  if (Gather(
          gather, values,
          [&dst, &ret, &type_owner](absl::Span<const npy_intp> size) {
            ret = PyArray_NewFromDescr(
                /*subtype=*/&PyArray_Type,
                /*type=*/
                reinterpret_cast<PyArray_Descr*>(type_owner.release().ptr()),
                /*nd=*/size.size(),
                /*dims=*/size.data(),
                /*strides=*/nullptr,
                /*data=*/nullptr,
                /*flags=*/0,
                /*obj=*/nullptr);
            if (ret == nullptr) {
              return false;
            }
            dst = static_cast<PyObject**>(
                PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)));
            return true;
          },
          [&dst, &existing, &last_value, &last_obj](absl::string_view value) {
            if (last_obj != nullptr && value == last_value) {
              Py_INCREF(last_obj);
              *dst++ = last_obj;
              return true;
            }
            auto [it, inserted] = existing.emplace(value, nullptr);
            if (inserted) {
              it->second =
                  PyUnicode_FromStringAndSize(value.data(), value.size());
              if (it->second == nullptr) {
                return false;
              }
            } else {
              Py_INCREF(it->second);
            }
            last_value = value;
            last_obj = it->second;
            *dst++ = last_obj;
            return true;
          })) {
    return ret;
  }
  Py_XDECREF(ret);
  return nullptr;
}

// Creates NumPy array with given dtype given specified converter.
// `converter` shall have the following signature:
// bool converter(const std::string& value, T* result);
// It must return whether conversion is successful and store conversion in
// result.
template <typename T, typename C>
inline PyObject* Convert(PyObject* gather, PyArray_Descr* type,
                         absl::Span<const std::string> values, C&& converter) {
  py::object ret;
  T* dst;
  std::optional<std::string> bad_value;
  if (Gather</*kReleaseGil=*/true>(
          gather, values,
          [&dst, &ret, type](absl::Span<const npy_intp> size) {
            // Construct uninitialised NumPy array of type T.
            ret = py::reinterpret_steal<py::object>(PyArray_NewFromDescr(
                /*subtype=*/&PyArray_Type,
                /*type=*/type,
                /*nd=*/size.size(),
                /*dims=*/size.data(),
                /*strides=*/nullptr,
                /*data=*/nullptr,
                /*flags=*/0,
                /*obj=*/nullptr));

            dst = static_cast<T*>(
                PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret.ptr())));
            return true;
          },
          [&dst, &converter, &bad_value](const std::string& value) {
            if (!converter(value, dst++)) {
              bad_value = value;
              return false;
            }
            return true;
          })) {
    return ret.release().ptr();
  }
  if (bad_value) {
    PyErr_SetString(PyExc_ValueError, bad_value->c_str());
  }
  return nullptr;
}

PyObject* CifDictGetArray(const CifDict& self, absl::string_view key,
                          PyObject* dtype, PyObject* gather) {
  import_array();
  PyArray_Descr* type = nullptr;
  if (dtype == Py_None) {
    type = PyArray_DescrFromType(NPY_OBJECT);
  } else if (PyArray_DescrConverter(dtype, &type) == NPY_FAIL || !type) {
    PyErr_Format(PyExc_TypeError, "Invalid dtype %R", dtype);
    Py_XDECREF(type);
    return nullptr;
  }
  auto entry = self.dict()->find(key);
  if (entry == self.dict()->end()) {
    Py_DECREF(type);
    PyErr_SetObject(PyExc_KeyError,
                    PyUnicode_FromStringAndSize(key.data(), key.size()));
    return nullptr;
  }

  auto int_convert = [](absl::string_view str, auto* value) {
    return absl::SimpleAtoi(str, value);
  };

  auto int_convert_bounded = [](absl::string_view str, auto* value) {
    int64_t v;
    if (absl::SimpleAtoi(str, &v)) {
      using limits =
          std::numeric_limits<std::remove_reference_t<decltype(*value)>>;
      if (limits::min() <= v && v <= limits::max()) {
        *value = v;
        return true;
      }
    }
    return false;
  };

  absl::Span<const std::string> values = entry->second;

  switch (type->type_num) {
    case NPY_DOUBLE:
      return Convert<double>(
          gather, type, values, [](absl::string_view str, double* value) {
            if (str == ".") {
              *value = std::numeric_limits<double>::quiet_NaN();
              return true;
            }
            return absl::SimpleAtod(str, value);
          });
    case NPY_FLOAT:
      return Convert<float>(
          gather, type, values, [](absl::string_view str, float* value) {
            if (str == ".") {
              *value = std::numeric_limits<float>::quiet_NaN();
              return true;
            }
            return absl::SimpleAtof(str, value);
          });
    case NPY_INT8:
      return Convert<int8_t>(gather, type, values, int_convert_bounded);
    case NPY_INT16:
      return Convert<int16_t>(gather, type, values, int_convert_bounded);
    case NPY_INT32:
      return Convert<int32_t>(gather, type, values, int_convert);
    case NPY_INT64:
      return Convert<int64_t>(gather, type, values, int_convert);
    case NPY_UINT8:
      return Convert<uint8_t>(gather, type, values, int_convert_bounded);
    case NPY_UINT16:
      return Convert<uint16_t>(gather, type, values, int_convert_bounded);
    case NPY_UINT32:
      return Convert<uint32_t>(gather, type, values, int_convert);
    case NPY_UINT64:
      return Convert<uint64_t>(gather, type, values, int_convert);
    case NPY_BOOL:
      return Convert<bool>(gather, type, values,
                           [](absl::string_view str, bool* value) {
                             if (str == "n" || str == "no") {
                               *value = false;
                               return true;
                             }
                             if (str == "y" || str == "yes") {
                               *value = true;
                               return true;
                             }
                             return false;
                           });
    case NPY_OBJECT:
      return ConvertStrings(gather, type, values);
    default: {
      PyErr_Format(PyExc_TypeError, "Unsupported dtype %R", dtype);
      Py_XDECREF(type);
      return nullptr;
    }
  }
}

}  // namespace

void RegisterModuleCifDict(pybind11::module m) {
  using Value = std::vector<std::string>;
  static absl::NoDestructor<std::vector<std::string>> empty_values;

  m.def(
      "from_string",
      [](absl::string_view s) {
        absl::StatusOr<CifDict> dict = CifDict::FromString(s);
        if (!dict.ok()) {
          throw py::value_error(dict.status().ToString());
        }
        return *dict;
      },
      py::call_guard<py::gil_scoped_release>());

  m.def(
      "tokenize",
      [](absl::string_view cif_string) {
        absl::StatusOr<std::vector<std::string>> tokens = Tokenize(cif_string);
        if (!tokens.ok()) {
          throw py::value_error(tokens.status().ToString());
        }
        return *std::move(tokens);
      },
      py::arg("cif_string"), py::call_guard<py::gil_scoped_release>());

  m.def("split_line", [](absl::string_view line) {
    absl::StatusOr<std::vector<absl::string_view>> tokens = SplitLine(line);
    if (!tokens.ok()) {
      throw py::value_error(tokens.status().ToString());
    }
    return *std::move(tokens);
  });

  m.def(
      "parse_multi_data_cif",
      [](absl::string_view cif_string) {
        absl::StatusOr<absl::flat_hash_map<std::string, CifDict>> result;
        {
          py::gil_scoped_release gil_release;
          result = ParseMultiDataCifDict(cif_string);
        }
        if (!result.ok()) {
          throw py::value_error(result.status().ToString());
        }
        py::dict dict;
        for (auto& [key, value] : *result) {
          dict[py::cast(key)] = py::cast(value);
        }
        return dict;
      },
      py::arg("cif_string"));

  auto cif_dict =
      py::class_<CifDict>(m, "CifDict")
          .def(py::init<>([](py::dict dict) {
                 CifDict::Dict result;
                 for (const auto& [key, value] : dict) {
                   result.emplace(py::cast<absl::string_view>(key),
                                  py::cast<std::vector<std::string>>(value));
                 }
                 return CifDict(std::move(result));
               }),
               "Initialise with a map")
          .def("copy_and_update",
               [](const CifDict& self, py::dict dict) {
                 CifDict::Dict result;
                 for (const auto& [key, value] : dict) {
                   result.emplace(py::cast<absl::string_view>(key),
                                  py::cast<std::vector<std::string>>(value));
                 }
                 {
                   py::gil_scoped_release gil_release;
                   return self.CopyAndUpdate(std::move(result));
                 }
               })
          .def(
              "__str__",
              [](const CifDict& self) {
                absl::StatusOr<std::string> result = self.ToString();
                if (!result.ok()) {
                  throw py::value_error(result.status().ToString());
                }
                return *result;
              },
              "Serialize to a string", py::call_guard<py::gil_scoped_release>())
          .def(
              "to_string",
              [](const CifDict& self) {
                absl::StatusOr<std::string> result = self.ToString();
                if (!result.ok()) {
                  throw py::value_error(result.status().ToString());
                }
                return *result;
              },
              "Serialize to a string", py::call_guard<py::gil_scoped_release>())
          .def(
              "to_dict",
              [](const CifDict& self) {
                py::dict result;
                for (const auto& [key, value] : *self.dict()) {
                  result[py::cast(key)] = py::cast(value);
                }
                return result;
              },
              "Returns the CIF data as a Python dict[str, list[str]].")
          .def("value_length", &CifDict::ValueLength, py::arg("key"),
               "Num elements in value")
          .def("__len__",
               [](const CifDict& self) { return self.dict()->size(); })
          .def(
              "__bool__",
              [](const CifDict& self) { return !self.dict()->empty(); },
              "Check whether the map is nonempty")
          .def(
              "__contains__",
              [](const CifDict& self, absl::string_view k) {
                return self.dict()->find(k) != self.dict()->end();
              },
              py::arg("key"), py::call_guard<py::gil_scoped_release>())
          .def("get_data_name", &CifDict::GetDataName)
          .def(
              "get",
              [](const CifDict& self, absl::string_view k,
                 py::object default_value) -> py::object {
                auto it = self.dict()->find(k);
                if (it == self.dict()->end()) return default_value;
                py::list result(it->second.size());
                size_t index = 0;
                for (const std::string& v : it->second) {
                  result[index++] = py::cast(v);
                }
                return result;
              },
              py::arg("key"), py::arg("default_value") = py::none())
          .def(
              "get_array",
              [](const CifDict& self, absl::string_view key, py::handle dtype,
                 py::handle gather) -> py::object {
                PyObject* obj =
                    CifDictGetArray(self, key, dtype.ptr(), gather.ptr());
                if (obj == nullptr) {
                  throw py::error_already_set();
                }
                return py::reinterpret_steal<py::object>(obj);
              },
              py::arg("key"), py::arg("dtype") = py::none(),
              py::arg("gather") = py::none())
          .def(
              "__getitem__",
              [](const CifDict& self, absl::string_view k) -> const Value& {
                auto it = self.dict()->find(k);
                if (it == self.dict()->end()) {
                  throw py::key_error(std::string(k).c_str());
                }
                return it->second;
              },
              py::arg("key"), py::call_guard<py::gil_scoped_release>())
          .def(
              "extract_loop_as_dict",
              [](const CifDict& self, absl::string_view prefix,
                 absl::string_view index) {
                absl::StatusOr<absl::flat_hash_map<
                    absl::string_view,
                    absl::flat_hash_map<absl::string_view, absl::string_view>>>
                    dict;
                {
                  py::gil_scoped_release gil_release;
                  dict = self.ExtractLoopAsDict(prefix, index);
                  if (!dict.ok()) {
                    throw py::value_error(dict.status().ToString());
                  }
                }
                py::dict key_value_dict;
                for (const auto& [key, value] : *dict) {
                  py::dict value_dict;
                  for (const auto& [key2, value2] : value) {
                    value_dict[py::cast(key2)] = py::cast(value2);
                  }
                  key_value_dict[py::cast(key)] = std::move(value_dict);
                }
                return key_value_dict;
              },
              py::arg("prefix"), py::arg("index"))
          .def(
              "extract_loop_as_list",
              [](const CifDict& self, absl::string_view prefix) {
                absl::StatusOr<std::vector<
                    absl::flat_hash_map<absl::string_view, absl::string_view>>>
                    list_dict;
                {
                  py::gil_scoped_release gil_release;
                  list_dict = self.ExtractLoopAsList(prefix);
                  if (!list_dict.ok()) {
                    throw py::value_error(list_dict.status().ToString());
                  }
                }
                py::list list_obj(list_dict->size());
                size_t index = 0;
                for (const auto& value : *list_dict) {
                  py::dict value_dict;
                  for (const auto& [key, value] : value) {
                    value_dict[py::cast(key)] = py::cast(value);
                  }
                  list_obj[index++] = std::move(value_dict);
                }
                return list_obj;
              },
              py::arg("prefix"))
          .def(py::pickle(
              [](const CifDict& self) {  // __getstate__.
                py::tuple result_tuple(1);
                py::dict result;
                for (const auto& [key, value] : *self.dict()) {
                  result[py::cast(key)] = py::cast(value);
                }
                result_tuple[0] = std::move(result);
                return result_tuple;
              },
              [](py::tuple t) {  // __setstate__.
                py::dict dict = t[0].cast<py::dict>();
                CifDict::Dict result;
                for (const auto& [key, value] : dict) {
                  result.emplace(py::cast<absl::string_view>(key),
                                 py::cast<std::vector<std::string>>(value));
                }
                return CifDict(std::move(result));
              }));

  // Item, value, and key views
  struct KeyView {
    CifDict map;
  };

  struct ValueView {
    CifDict map;
  };
  struct ItemView {
    CifDict map;
  };

  py::class_<ItemView>(cif_dict, "ItemView")
      .def("__len__", [](const ItemView& v) { return v.map.dict()->size(); })
      .def(
          "__iter__",
          [](const ItemView& v) {
            return py::make_iterator(v.map.dict()->begin(),
                                     v.map.dict()->end());
          },
          py::keep_alive<0, 1>());

  py::class_<KeyView>(cif_dict, "KeyView")
      .def(
          "__contains__",
          [](const KeyView& v, absl::string_view k) {
            return v.map.dict()->find(k) != v.map.dict()->end();
          },
          py::call_guard<py::gil_scoped_release>())
      .def("__contains__", [](const KeyView&, py::handle) { return false; })
      .def("__len__", [](const KeyView& v) { return v.map.dict()->size(); })
      .def(
          "__iter__",
          [](const KeyView& v) {
            return py::make_key_iterator(v.map.dict()->begin(),
                                         v.map.dict()->end());
          },
          py::keep_alive<0, 1>());

  py::class_<ValueView>(cif_dict, "ValueView")
      .def("__len__", [](const ValueView& v) { return v.map.dict()->size(); })
      .def(
          "__iter__",
          [](const ValueView& v) {
            return py::make_value_iterator(v.map.dict()->begin(),
                                           v.map.dict()->end());
          },
          py::keep_alive<0, 1>());

  cif_dict
      .def(
          "__iter__",
          [](const CifDict& self) {
            return py::make_key_iterator(self.dict()->begin(),
                                         self.dict()->end());
          },
          py::keep_alive<0, 1>())
      .def(
          "keys", [](const CifDict& self) { return KeyView{self}; },
          "Returns an iterable view of the map's keys.")
      .def(
          "values", [](const CifDict& self) { return ValueView{self}; },
          "Returns an iterable view of the map's values.")
      .def(
          "items", [](const CifDict& self) { return ItemView{self}; },
          "Returns an iterable view of the map's items.");
}

}  // namespace alphafold3
