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

/// Project inlcudes
#include "queso/python/bindings/define_python.hpp"
#include "queso/python/bindings/add_dictionary_to_python.h"
#include "queso/python/bindings/dictionary_binder_helper.hpp"
// To export
#include "queso/includes/dictionary_factory.hpp"
#include "queso/containers/dictionary.hpp"

// Note: PYBIND11_MAKE_OPAQUE can not be captured within namespace
using DictionaryType = queso::Dictionary<queso::key::MainValuesTypeTag>;
using DictionaryPtrType = queso::Unique<DictionaryType>;
using DictionaryVectorPtrType = std::vector<DictionaryPtrType>;
PYBIND11_MAKE_OPAQUE(DictionaryVectorPtrType);

namespace queso {
namespace Python {

using MainDictionaryHolderType = UniqueHolder<MainDictionaryType>;

namespace py = pybind11;

void AddDictionaryToPython(pybind11::module& m) {

    // Export DictionaryHolder
    py::class_<MainDictionaryHolderType>(m,"DictionaryHolder")
        .def("GetObject", &MainDictionaryHolderType::GetObject, py::return_value_policy::reference_internal)
    ;

    // Export Dictionary
    py::class_<DictionaryType, DictionaryPtrType> dictionary_base_type_binder(m,"Dictionary");

    // Export methods
    DictionaryBinderHelper<DictionaryType>(dictionary_base_type_binder);


    // Export static methods
    dictionary_base_type_binder.def_static("Create", [](const std::string& rName) -> MainDictionaryHolderType {
        return MainDictionaryHolderType(DictionaryFactory<DictionaryType::KeySetValuesTypeTag>::Create(rName)); },
            py::return_value_policy::move );

    // Export vector
    py::class_<DictionaryVectorPtrType>(m, "DictionaryList")
        .def("__getitem__", [](DictionaryVectorPtrType &self, IndexType i)
            { return &(*self[i]); }, py::return_value_policy::reference_internal)
        .def("__len__", [](const DictionaryVectorPtrType &self) { return self.size(); })
        .def("__iter__", [](DictionaryVectorPtrType &self) {
            return py::make_iterator( dereference_iterator(self.begin()), dereference_iterator(self.end()) );
        }, py::keep_alive<0, 1>() )
        .def("append", [](DictionaryVectorPtrType& self, MainDictionaryHolderType& rDictionaryHolder){
            self.push_back(std::move(rDictionaryHolder.Release()));
        } )
    ;

} // End AddDictionaryToPython

} // End namespace Python
} // End namespace queso

