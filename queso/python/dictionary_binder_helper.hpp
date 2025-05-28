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

#ifndef DICTIONARY_BINDER_HELPER_INCLUDE_HPP
#define DICTIONARY_BINDER_HELPER_INCLUDE_HPP

/// Project inlcudes
#include "queso/python/define_python.hpp"
#include "queso/containers/dictionary.hpp"

namespace queso {
namespace Python {

namespace py = pybind11;

template<class TDictType>
void DictionaryBinderHelper(py::class_<TDictType, Unique<TDictType>>& binder) {

    using StringAccess = DictionaryStringAccess<TDictType>;
    /// Print
    binder.def("__str__", PrintObject<TDictType>);

    /// Get sub dictionaries
    binder.def("__getitem__", [](TDictType& rDictionary, const std::string& rKeyName) {
        return &StringAccess::GetSubDictionary(rDictionary, rKeyName); },
            py::return_value_policy::reference_internal );

    /// Get list
    binder.def("GetList", [](TDictType& rDictionary, const std::string& rKeyName){
        return &StringAccess::GetList(rDictionary, rKeyName);},
            py::return_value_policy::reference_internal );

    /// IsSet
    binder.def("IsSet", [](const TDictType& rDictionary, const std::string& rKeyName){
        return StringAccess::IsSet(rDictionary, rKeyName); } );

    /// SetValue
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, const PointType& rValue){
        StringAccess::SetValue(rDictionary, rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, const Vector3i& rValue){
        StringAccess::SetValue(rDictionary, rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, bool rValue){
        StringAccess::SetValue(rDictionary, rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, double rValue){
        StringAccess::SetValue(rDictionary, rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, IndexType rValue){
        StringAccess::SetValue(rDictionary, rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, const std::string& rValue){
        StringAccess::SetValue(rDictionary, rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, IntegrationMethodType rValue){
        StringAccess::SetValue(rDictionary, rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, const GridTypeType& rValue){
        StringAccess::SetValue(rDictionary, rKeyName, rValue); });

    /// GetValue
    binder.def("GetDoubleVector", [](const TDictType& rDictionary, const std::string& rKeyName) -> PointType {
        return StringAccess::template GetValue<PointType>(rDictionary, rKeyName); });
    binder.def("GetIntVector", [](const TDictType& rDictionary, const std::string& rKeyName ) -> Vector3i{
        return StringAccess::template GetValue<Vector3i>(rDictionary, rKeyName); });
    binder.def("GetDouble", [](const TDictType& rDictionary, const std::string& rKeyName ) -> double {
        return StringAccess::template GetValue<double>(rDictionary, rKeyName); });
    binder.def("GetBool", [](const TDictType& rDictionary, const std::string& rKeyName ) -> bool {
        return StringAccess::template GetValue<bool>(rDictionary, rKeyName); });
    binder.def("GetInt", [](const TDictType& rDictionary, const std::string& rKeyName ) -> IndexType {
        return StringAccess::template GetValue<IndexType>(rDictionary, rKeyName); });
    binder.def("GetString", [](const TDictType& rDictionary, const std::string& rKeyName ) -> std::string {
        return StringAccess::template GetValue<std::string>(rDictionary, rKeyName); });
    binder.def("GetIntegrationMethod", [](const TDictType& rDictionary, const std::string& rKeyName ) -> IntegrationMethodType {
        return StringAccess::template GetValue<IntegrationMethodType>(rDictionary, rKeyName); });
    binder.def("GetGridType", [](const TDictType& rDictionary, const std::string& rKeyName ) -> GridTypeType {
        return StringAccess::template GetValue<GridTypeType>(rDictionary, rKeyName); });
}

} // End Python
} // End queso

#endif // DICTIONARY_BINDER_HELPER_INCLUDE_HPP