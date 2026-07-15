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

#include <utility>

#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "alphafold3/parsers/cpp/fasta_iterator_lib.h"
#include "pybind11/attr.h"
#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"
#include "pybind11/stl.h"

namespace alphafold3 {
namespace {

namespace py = pybind11;

template <typename T>
T ValueOrThrowValueError(absl::StatusOr<T> value) {
  if (!value.ok()) throw py::value_error(value.status().ToString());
  return *std::move(value);
}

constexpr char kFastaFileIteratorDoc[] = R"(
Lazy FASTA parser for memory efficient FASTA parsing from a path.)";

constexpr char kFastaStringIteratorDoc[] = R"(
Lazy FASTA parser for memory efficient FASTA parsing from a string.

WARNING: The object backing the fasta_string string_view must not be
deleted while the FastaStringIterator is alive. E.g. this will break:

```
# Make sure the fasta_string is not interned.
fasta_string = '\n'.join(['>d\nS' for _ in range(10)])
iterator = fasta_iterator.FastaStringIterator(fasta_string)
del fasta_string
iterator.next()  # Heap use-after-free.
```
)";

constexpr char kParseFastaDoc[] = R"(
Parses a FASTA string and returns a list of amino-acid sequences.

Args:
  fasta_string: The contents of a FASTA file.

Returns:
  List of sequences in the FASTA file. Descriptions are ignored.
)";

constexpr char kParseFastaIncludeDescriptionsDoc[] = R"(
Parses a FASTA string, returns amino-acid sequences with descriptions.

Args:
  fasta_string: The contents of a FASTA file.

Returns:
  A tuple with two lists (sequences, descriptions):
  * A list of sequences.
  * A list of sequence descriptions taken from the comment lines. In the
    same order as the sequences.
)";

class PythonFastaStringIterator : public FastaStringIterator {
 public:
  explicit PythonFastaStringIterator(py::object fasta_string)
      : FastaStringIterator(py::cast<absl::string_view>(fasta_string)),
        fasta_string_(std::move(fasta_string)) {}

 private:
  py::object fasta_string_;
};

}  // namespace

void RegisterModuleFastaIterator(pybind11::module m) {
  py::class_<FastaFileIterator>(m, "FastaFileIterator", kFastaFileIteratorDoc)
      .def(py::init<absl::string_view>(), py::arg("fasta_path"))
      .def("__iter__",
           [](FastaFileIterator& iterator) -> FastaFileIterator& {
             return iterator;
           })
      .def(
          "__next__",
          [](FastaFileIterator& iterator) {
            if (iterator.HasNext()) {
              return ValueOrThrowValueError(iterator.Next());
            } else {
              throw py::stop_iteration();
            }
          },
          py::call_guard<py::gil_scoped_release>());

  py::class_<PythonFastaStringIterator>(m, "FastaStringIterator",
                                        kFastaStringIteratorDoc)
      .def(py::init<py::object>(), py::arg("fasta_string"))
      .def("__iter__",
           [](PythonFastaStringIterator& iterator)
               -> PythonFastaStringIterator& { return iterator; })
      .def(
          "__next__",
          [](PythonFastaStringIterator& iterator) {
            if (iterator.HasNext()) {
              return ValueOrThrowValueError(iterator.Next());
            } else {
              throw py::stop_iteration();
            }
          },
          py::call_guard<py::gil_scoped_release>());

  m.def("parse_fasta", &ParseFasta, py::arg("fasta_string"),
        py::call_guard<py::gil_scoped_release>(), py::doc(kParseFastaDoc + 1));
  m.def("parse_fasta_include_descriptions", &ParseFastaIncludeDescriptions,
        py::arg("fasta_string"), py::call_guard<py::gil_scoped_release>(),
        py::doc(kParseFastaIncludeDescriptionsDoc + 1));
}

}  // namespace alphafold3
