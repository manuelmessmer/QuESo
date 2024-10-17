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
void DictionaryBinderHelper(py::class_<TDictType>& binder) {

    /// Print
    binder.def("__str__", PrintObject<TDictType>);

    /// Get sub dictionaries
    binder.def("__getitem__", [](TDictType& rDictionary, const std::string& rKeyName){
        return &rDictionary[rKeyName]; }, py::return_value_policy::reference_internal );

    /// Get list
    binder.def("GetList", [](TDictType& rDictionary, const std::string& rKeyName){
        return &rDictionary.GetList(rKeyName); }, py::return_value_policy::reference_internal );

    /// IsSet
    binder.def("IsSet", [](const TDictType& rDictionary, const std::string& rKeyName){
        return rDictionary.IsSet(rKeyName); } );

    /// SetValue
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, const PointType& rValue){
        rDictionary.SetValue(rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, const Vector3i& rValue){
        rDictionary.SetValueWithAmbiguousType(rKeyName, rValue); }); // Allows cast to PointType, if possible.
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, bool rValue){
        rDictionary.SetValue(rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, double rValue){
        rDictionary.SetValue(rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, IndexType rValue){
        rDictionary.SetValueWithAmbiguousType(rKeyName, rValue); }); // Allows cast to double, if possible.
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, const std::string& rValue){
        rDictionary.SetValue(rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, IntegrationMethodType rValue){
        rDictionary.SetValue(rKeyName, rValue); });
    binder.def("SetValue", [](TDictType& rDictionary, const std::string& rKeyName, const GridTypeType& rValue){
        rDictionary.SetValue(rKeyName, rValue); });

    /// GetValue
    binder.def("GetDoubleVector", [](TDictType& rDictionary, const std::string& rKeyName) -> PointType {
        return rDictionary.template GetValue<PointType>(rKeyName); });
    binder.def("GetIntVector", [](const TDictType& rDictionary, const std::string& rKeyName ) -> Vector3i{
        return rDictionary.template GetValue<Vector3i>(rKeyName); });
    binder.def("GetDouble", [](const TDictType& rDictionary, const std::string& rKeyName ) -> double {
        return rDictionary.template GetValue<double>(rKeyName); });
    binder.def("GetBool", [](const TDictType& rDictionary, const std::string& rKeyName ) -> bool {
        return rDictionary.template GetValue<bool>(rKeyName); });
    binder.def("GetInt", [](const TDictType& rDictionary, const std::string& rKeyName ) -> IndexType {
        return rDictionary.template GetValue<IndexType>(rKeyName); });
    binder.def("GetString", [](const TDictType& rDictionary, const std::string& rKeyName ) -> std::string {
        return rDictionary.template GetValue<std::string>(rKeyName); });
    binder.def("GetIntegrationMethod", [](const TDictType& rDictionary, const std::string& rKeyName ) -> IntegrationMethodType {
        return rDictionary.template GetValue<IntegrationMethodType>(rKeyName); });
    binder.def("GetGridType", [](const TDictType& rDictionary, const std::string& rKeyName ) -> GridTypeType {
        return rDictionary.template GetValue<GridTypeType>(rKeyName); });
}

} // End Python
} // End queso

#endif // DICTIONARY_BINDER_HELPER_INCLUDE_HPP