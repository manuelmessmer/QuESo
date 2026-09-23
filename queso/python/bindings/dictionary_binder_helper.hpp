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

#pragma once

// Project includes
#include "queso/containers/dictionary.hpp"
#include "queso/python/bindings/define_python.hpp"

namespace queso::python {

namespace py = pybind11;

template<class TDictType>
void DictionaryBinderHelper(py::class_<TDictType, Unique<TDictType>>& binder)
{
    using namespace pybind11::literals;
    using StringAccess = DictionaryStringAccess<TDictType>;
    binder.def("__str__", PrintObject<TDictType>, "Return a human-readable representation.");

    binder.def(
        "__getitem__",
        [](TDictType& rDictionary, const std::string& rKeyName) {
            return &StringAccess::GetSubDictionary(rDictionary, rKeyName);
        },
        py::return_value_policy::reference_internal,
        "Return a subdictionary by key."
    );

    binder.def(
        "get_list",
        [](TDictType& rDictionary, const std::string& rKeyName) {
            return &StringAccess::GetList(rDictionary, rKeyName);
        },
        "key"_a,
        py::return_value_policy::reference_internal,
        "Return a dictionary list by key."
    );

    binder.def(
        "is_set",
        [](const TDictType& rDictionary, const std::string& rKeyName) {
            return StringAccess::IsSet(rDictionary, rKeyName);
        },
        "key"_a,
        "Return whether a value is set."
    );

    binder.def(
        "set_value",
        [](TDictType& rDictionary, const std::string& rKeyName, const PointType& rValue) {
            StringAccess::SetValue(rDictionary, rKeyName, rValue);
        },
        "key"_a,
        "value"_a,
        "Set a three-component floating-point value."
    );
    binder.def(
        "set_value",
        [](TDictType& rDictionary, const std::string& rKeyName, const Vector3i& rValue) {
            StringAccess::SetValue(rDictionary, rKeyName, rValue);
        },
        "key"_a,
        "value"_a,
        "Set a three-component integer value."
    );
    binder.def(
        "set_value",
        [](TDictType& rDictionary, const std::string& rKeyName, bool rValue) {
            StringAccess::SetValue(rDictionary, rKeyName, rValue);
        },
        "key"_a,
        "value"_a,
        "Set a boolean value."
    );
    binder.def(
        "set_value",
        [](TDictType& rDictionary, const std::string& rKeyName, double rValue) {
            StringAccess::SetValue(rDictionary, rKeyName, rValue);
        },
        "key"_a,
        "value"_a,
        "Set a floating-point value."
    );
    binder.def(
        "set_value",
        [](TDictType& rDictionary, const std::string& rKeyName, IndexType rValue) {
            StringAccess::SetValue(rDictionary, rKeyName, rValue);
        },
        "key"_a,
        "value"_a,
        "Set an integer value."
    );
    binder.def(
        "set_value",
        [](TDictType& rDictionary, const std::string& rKeyName, const std::string& rValue) {
            StringAccess::SetValue(rDictionary, rKeyName, rValue);
        },
        "key"_a,
        "value"_a,
        "Set a string value."
    );
    binder.def(
        "set_value",
        [](TDictType& rDictionary, const std::string& rKeyName, IntegrationMethodType rValue) {
            StringAccess::SetValue(rDictionary, rKeyName, rValue);
        },
        "key"_a,
        "value"_a,
        "Set an integration method."
    );
    binder.def(
        "set_value",
        [](TDictType& rDictionary, const std::string& rKeyName, const GridTypeType& rValue) {
            StringAccess::SetValue(rDictionary, rKeyName, rValue);
        },
        "key"_a,
        "value"_a,
        "Set a background-grid type."
    );

    binder.def(
        "get_double_vector",
        [](const TDictType& rDictionary, const std::string& rKeyName) -> PointType {
            return StringAccess::template GetValue<PointType>(rDictionary, rKeyName);
        },
        "key"_a,
        "Return a three-component floating-point value."
    );
    binder.def(
        "get_int_vector",
        [](const TDictType& rDictionary, const std::string& rKeyName) -> Vector3i {
            return StringAccess::template GetValue<Vector3i>(rDictionary, rKeyName);
        },
        "key"_a,
        "Return a three-component integer value."
    );
    binder.def(
        "get_double",
        [](const TDictType& rDictionary, const std::string& rKeyName) -> double {
            return StringAccess::template GetValue<double>(rDictionary, rKeyName);
        },
        "key"_a,
        "Return a floating-point value."
    );
    binder.def(
        "get_bool",
        [](const TDictType& rDictionary, const std::string& rKeyName) -> bool {
            return StringAccess::template GetValue<bool>(rDictionary, rKeyName);
        },
        "key"_a,
        "Return a boolean value."
    );
    binder.def(
        "get_int",
        [](const TDictType& rDictionary, const std::string& rKeyName) -> IndexType {
            return StringAccess::template GetValue<IndexType>(rDictionary, rKeyName);
        },
        "key"_a,
        "Return an integer value."
    );
    binder.def(
        "get_string",
        [](const TDictType& rDictionary, const std::string& rKeyName) -> std::string {
            return StringAccess::template GetValue<std::string>(rDictionary, rKeyName);
        },
        "key"_a,
        "Return a string value."
    );
    binder.def(
        "get_integration_method",
        [](const TDictType& rDictionary, const std::string& rKeyName) -> IntegrationMethodType {
            return StringAccess::template GetValue<IntegrationMethodType>(rDictionary, rKeyName);
        },
        "key"_a,
        "Return an integration method."
    );
    binder.def(
        "get_grid_type",
        [](const TDictType& rDictionary, const std::string& rKeyName) -> GridTypeType {
            return StringAccess::template GetValue<GridTypeType>(rDictionary, rKeyName);
        },
        "key"_a,
        "Return a background-grid type."
    );
}

}  // namespace queso::python
