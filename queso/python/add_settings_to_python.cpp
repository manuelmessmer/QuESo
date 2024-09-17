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
#include "queso/python/define_python.hpp"
#include "queso/python/add_settings_to_python.h"
// To export
#include "queso/includes/settings.hpp"

namespace queso {
namespace Python {

namespace py = pybind11;

template< class TBinderType, typename TEnumType>
void DictionaryBinderHelper(TBinderType& binder) {

    binder.def("__getitem__", [](SettingsBaseType& rDictionary, TEnumType rEnum){
        return &rDictionary[rEnum]; }, py::return_value_policy::reference_internal );
    binder.def("__str__", PrintObject<SettingsBaseType>);
    binder.def("IsSet", [](SettingsBaseType& rDictionary, TEnumType rEnum){
        return rDictionary.IsSet(rEnum); } );

    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, const PointType& rValue){
        rDictionary.SetValue(rEnum, rValue); });
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, const Vector3i& rValue){
        rDictionary.SetValueWithAmbiguousType(rEnum, rValue); }); // Allows cast to PointType, if possible.
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, bool rValue){
        rDictionary.SetValue(rEnum, rValue); });
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, double rValue){
        rDictionary.SetValue(rEnum, rValue); });
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, IndexType rValue){
        rDictionary.SetValueWithAmbiguousType(rEnum, rValue); }); // Allows cast to double, if possible.
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, const std::string& rValue){
        rDictionary.SetValue(rEnum, rValue); });
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, IntegrationMethodType rValue){
        rDictionary.SetValue(rEnum, rValue); });
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, const BackgroundGridTypeType& rValue){
        rDictionary.SetValue(rEnum, rValue); });

    binder.def("GetDoubleVector", [](SettingsBaseType& rDictionary, TEnumType rEnum){
        return rDictionary.GetValue<PointType>(rEnum); });
    binder.def("GetIntVector", [](SettingsBaseType& rDictionary, TEnumType rEnum ){
        return rDictionary.GetValue<Vector3i>(rEnum); });
    binder.def("GetDouble", [](SettingsBaseType& rDictionary, TEnumType rEnum ){
        return rDictionary.GetValue<double>(rEnum); });
    binder.def("GetBool", [](SettingsBaseType& rDictionary, TEnumType rEnum ){
        return rDictionary.GetValue<bool>(rEnum); });
    binder.def("GetInt", [](SettingsBaseType& rDictionary, TEnumType rEnum ){
        return rDictionary.GetValue<IndexType>(rEnum); });
    binder.def("GetString", [](SettingsBaseType& rDictionary, TEnumType rEnum ){
        return rDictionary.GetValue<std::string>(rEnum); });
    binder.def("GetIntegrationMethod", [](SettingsBaseType& rDictionary, TEnumType rEnum ){
        return rDictionary.GetValue<IntegrationMethodType>(rEnum); });
    binder.def("GetBackgroundGridType", [](SettingsBaseType& rDictionary, TEnumType rEnum ){
        return rDictionary.GetValue<BackgroundGridTypeType>(rEnum); });
}

void AddSettingsToPython(pybind11::module& m) {

    /// Export enum MainSettings
    py::enum_<MainSettings> (m, "MainSettings")
        .value("general_settings", MainSettings::general_settings)
        .value("background_grid_settings", MainSettings::background_grid_settings)
        .value("trimmed_quadrature_rule_settings", MainSettings::trimmed_quadrature_rule_settings)
        .value("non_trimmed_quadrature_rule_settings", MainSettings::non_trimmed_quadrature_rule_settings)
        .value("conditions_settings_list", MainSettings::conditions_settings_list)
    ;

    /// Export enum GeneralSettings
    py::enum_<GeneralSettings> (m, "GeneralSettings")
        .value("input_filename", GeneralSettings::input_filename)
        .value("output_directory_name", GeneralSettings::output_directory_name)
        .value("echo_level", GeneralSettings::echo_level)
        .value("write_output_to_file", GeneralSettings::write_output_to_file)
    ;

    /// Export enum BackgroundGridSettings
    py::enum_<BackgroundGridSettings> (m, "BackgroundGridSettings")
        .value("grid_type", BackgroundGridSettings::grid_type)
        .value("lower_bound_xyz", BackgroundGridSettings::lower_bound_xyz)
        .value("upper_bound_xyz", BackgroundGridSettings::upper_bound_xyz)
        .value("lower_bound_uvw", BackgroundGridSettings::lower_bound_uvw)
        .value("upper_bound_uvw", BackgroundGridSettings::upper_bound_uvw)
        .value("polynomial_order", BackgroundGridSettings::polynomial_order)
        .value("number_of_elements", BackgroundGridSettings::number_of_elements)
    ;

    /// Export enum TrimmedQuadratureRuleSettings
    py::enum_<TrimmedQuadratureRuleSettings> (m, "TrimmedQuadratureRuleSettings")
        .value("moment_fitting_residual", TrimmedQuadratureRuleSettings::moment_fitting_residual)
        .value("min_element_volume_ratio", TrimmedQuadratureRuleSettings::min_element_volume_ratio)
        .value("min_num_boundary_triangles", TrimmedQuadratureRuleSettings::min_num_boundary_triangles)
        .value("neglect_elements_if_stl_is_flawed", TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed)
    ;

    /// Export enum NonTrimmedQuadratureRuleSettings
    py::enum_<NonTrimmedQuadratureRuleSettings> (m, "NonTrimmedQuadratureRuleSettings")
        .value("integration_method", NonTrimmedQuadratureRuleSettings::integration_method)
    ;

    /// Export enum ConditionSettings
    py::enum_<ConditionSettings> (m, "ConditionSettings")
        .value("condition_id", ConditionSettings::condition_id)
        .value("condition_type", ConditionSettings::condition_type)
        .value("input_filename", ConditionSettings::input_filename)
        .value("modulus", ConditionSettings::modulus)
        .value("direction", ConditionSettings::direction)
        .value("value", ConditionSettings::value)
        .value("penalty_factor", ConditionSettings::penalty_factor)
    ;

    /// Export Dictionary
    typedef py::class_<SettingsBaseType> SettingsBaseTypeBinderType;
    py::class_<SettingsBaseType> settings_base_type_binder(m,"Dictionary");

    DictionaryBinderHelper<SettingsBaseTypeBinderType, MainSettings>(settings_base_type_binder);
    DictionaryBinderHelper<SettingsBaseTypeBinderType, GeneralSettings>(settings_base_type_binder);
    DictionaryBinderHelper<SettingsBaseTypeBinderType, BackgroundGridSettings>(settings_base_type_binder);
    DictionaryBinderHelper<SettingsBaseTypeBinderType, TrimmedQuadratureRuleSettings>(settings_base_type_binder);
    DictionaryBinderHelper<SettingsBaseTypeBinderType, NonTrimmedQuadratureRuleSettings>(settings_base_type_binder);
    DictionaryBinderHelper<SettingsBaseTypeBinderType, ConditionSettings>(settings_base_type_binder);

    settings_base_type_binder.def("NumberOfSubDictionaries", &SettingsBaseType::NumberOfSubDictionaries);
    settings_base_type_binder.def("__getitem__", [](SettingsBaseType& rDictionary, IndexType Index){
        return &rDictionary[Index]; }, py::return_value_policy::reference_internal );
    settings_base_type_binder.def("__iter__",    [](SettingsBaseType& rDictionary){
        return py::make_iterator(rDictionary.begin_sub_dicts(), rDictionary.end_sub_dicts());},  py::keep_alive<0,1>());
    settings_base_type_binder.def("__len__",    [](SettingsBaseType& rDictionary){return rDictionary.NumberOfSubDictionaries(); });

    /// Export Settings
    py::class_<Settings, SettingsBaseType>(m,"Settings")
        .def(py::init<>())
        .def("CreateNewConditionSettings", &Settings::CreateNewConditionSettings, py::return_value_policy::reference_internal)
        .def("Check", &Settings::Check)
    ;

} // End AddSettingsToPython

} // End namespace Python
} // End namespace queso

