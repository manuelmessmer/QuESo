//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

// Project includes
#include "queso/python/bindings/add_dictionary_to_python.h"
#include "queso/containers/dictionary.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/python/bindings/define_python.hpp"
#include "queso/python/bindings/dictionary_binder_helper.hpp"

// Note: PYBIND11_MAKE_OPAQUE can not be captured within namespace
using DictionaryType = queso::Dictionary<queso::key::MainValuesTypeTag>;
using DictionaryPtrType = queso::Unique<DictionaryType>;
using DictionaryVectorPtrType = std::vector<DictionaryPtrType>;
PYBIND11_MAKE_OPAQUE(DictionaryVectorPtrType);

namespace queso::python {

using MainDictionaryHolderType = UniqueHolder<MainDictionaryType>;

namespace py = pybind11;
using namespace pybind11::literals;

void AddDictionaryToPython(pybind11::module& m)
{
    py::class_<MainDictionaryHolderType>(m, "DictionaryHolder", "Temporary owner used to transfer a dictionary.")
        .def_property_readonly(
            "dictionary",
            &MainDictionaryHolderType::GetObject,
            py::return_value_policy::reference_internal,
            "Dictionary owned by this holder."
        );

    py::class_<DictionaryType, DictionaryPtrType> DictionaryBinder(
        m, "Dictionary", "Hierarchical settings dictionary."
    );

    DictionaryBinderHelper<DictionaryType>(DictionaryBinder);
    DictionaryBinder.def_static(
        "create",
        [](const std::string& rName) -> MainDictionaryHolderType {
            return MainDictionaryHolderType(DictionaryFactory<DictionaryType::KeySetValuesTypeTag>::Create(rName));
        },
        "name"_a,
        py::return_value_policy::move,
        "Create a dictionary with the named schema."
    );

    py::class_<DictionaryVectorPtrType>(
        m, "DictionaryList", "Mutable list of dictionaries owned by a parent dictionary."
    )
        .def(
            "__getitem__",
            [](DictionaryVectorPtrType& self, py::ssize_t Index) {
                const auto size = static_cast<py::ssize_t>(self.size());
                if (Index < 0) { Index += size; }
                if (Index < 0 || Index >= size) { throw py::index_error(); }
                return self[static_cast<IndexType>(Index)].get();
            },
            py::return_value_policy::reference_internal,
            "Return a dictionary by index."
        )
        .def(
            "__len__", [](const DictionaryVectorPtrType& self) { return self.size(); }, "Number of dictionaries."
        )
        .def(
            "__iter__",
            [](DictionaryVectorPtrType& self) {
                return py::make_iterator(dereference_iterator(self.begin()), dereference_iterator(self.end()));
            },
            py::keep_alive<0, 1>()
        )
        .def(
            "append",
            [](DictionaryVectorPtrType& self, MainDictionaryHolderType& rDictionaryHolder) {
                if (!rDictionaryHolder.HasObject()) {
                    throw py::value_error("DictionaryHolder no longer owns a dictionary.");
                }
                self.push_back(std::move(rDictionaryHolder.Release()));
            },
            "dictionary_holder"_a,
            "Move a dictionary from ``dictionary_holder`` into this list. The holder is empty afterward."
        );
}

}  // namespace queso::python
