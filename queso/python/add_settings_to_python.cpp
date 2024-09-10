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

/// External includes
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
/// Project inlcudes
#include "queso/python/add_settings_to_python.h"
#include "queso/includes/settings.hpp"

namespace queso {
namespace Python {

namespace py = pybind11;

template <typename TType>
TType pyStringToEnum(const py::enum_<TType>& enm, const std::string& value) {
    const auto values = enm.attr("__members__").template cast<py::dict>();
    const auto strVal = py::str(value);
    if (values.contains(strVal)) {
        return TType(values[strVal].template cast<TType>());
    }
    QuESo_ERROR << "Invalid string value " << value << " for enum " << std::string(typeid(TType).name());
}
//general_settings, mesh_settings, trimmed_quadrature_rule_settings, non_trimmed_quadrature_rule_settings, conditions_settings_list
static std::map<std::string, Main> main_map = {
    {std::string("general_settings"), Main::general_settings},
    {std::string("mesh_settings"), Main::mesh_settings},
    {std::string("trimmed_quadrature_rule_settings"), Main::trimmed_quadrature_rule_settings},
    {std::string("non_trimmed_quadrature_rule_settings"), Main::non_trimmed_quadrature_rule_settings},
    {std::string("conditions_settings_list"), Main::conditions_settings_list}
};

template< class TBinderType, typename TEnumType>
void DictionaryBinderHelper(TBinderType& binder) {

    binder.def("__getitem__", [](SettingsBaseType& rDictionary, TEnumType rEnum){
        return &rDictionary[rEnum]; }, py::return_value_policy::reference_internal );
    binder.def("IsSet", [](SettingsBaseType& rDictionary, TEnumType rEnum){
        return rDictionary.IsSet(rEnum); } );

    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, const PointType& rValue){
        rDictionary.SetValue(rEnum, rValue); });
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, const Vector3i& rValue){
        rDictionary.SetValue(rEnum, rValue); });
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, bool rValue){
        rDictionary.SetValue(rEnum, rValue); });
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, unsigned long rValue){
        rDictionary.SetValue(rEnum, rValue); });
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, IntegrationMethodType rValue){
        rDictionary.SetValue(rEnum, rValue); });
    binder.def("SetValue", [](SettingsBaseType& rDictionary, TEnumType rEnum, const std::string& rValue){
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
}



void AddSettingsToPython(pybind11::module& m) {

    /// Export enum IntegrationMethod
    py::enum_<Main> enm_main(m, "Main");
    enm_main
        .value("general_settings", Main::general_settings)
        .value("mesh_settings", Main::mesh_settings)
        .value("trimmed_quadrature_rule_settings", Main::trimmed_quadrature_rule_settings)
        .value("non_trimmed_quadrature_rule_settings", Main::non_trimmed_quadrature_rule_settings)
        .value("conditions_settings_list", Main::conditions_settings_list)
    ;
    py::implicitly_convertible<std::string, Main>();

    py::enum_<GeneralSettings> enm_general_settings(m, "GeneralSettings");
    enm_general_settings
        .value("input_filename", GeneralSettings::input_filename)
        .value("output_directory_name", GeneralSettings::output_directory_name)
        .value("echo_level", GeneralSettings::echo_level)
        .value("embedding_flag", GeneralSettings::embedding_flag)
    ;
    py::implicitly_convertible<std::string, GeneralSettings>();

    py::enum_<MeshSettings> enm_mesh_settings(m, "MeshSettings");
    enm_mesh_settings
        .value("lower_bound_xyz", MeshSettings::lower_bound_xyz)
        .value("upper_bound_xyz", MeshSettings::upper_bound_xyz)
        .value("lower_bound_uvw", MeshSettings::lower_bound_uvw)
        .value("upper_bound_uvw", MeshSettings::upper_bound_uvw)
        .value("polynomial_order", MeshSettings::polynomial_order)
        .value("number_of_elements", MeshSettings::number_of_elements)
    ;
    py::implicitly_convertible<std::string, MeshSettings>();

    py::enum_<TrimmedQuadratureRuleSettings> enm_trimmed_quad_rule_settings(m, "TrimmedQuadratureRuleSettings");
    enm_trimmed_quad_rule_settings
        .value("moment_fitting_residual", TrimmedQuadratureRuleSettings::moment_fitting_residual)
        .value("min_element_volume_ratio", TrimmedQuadratureRuleSettings::min_element_volume_ratio)
        .value("min_num_boundary_triangles", TrimmedQuadratureRuleSettings::min_num_boundary_triangles)
        .value("use_customized_trimmed_points", TrimmedQuadratureRuleSettings::use_customized_trimmed_points)
    ;
    py::implicitly_convertible<std::string, TrimmedQuadratureRuleSettings>();

    py::enum_<NonTrimmedQuadratureRuleSettings> enm_non_trimmed_quad_rule_settings(m, "NonTrimmedQuadratureRuleSettings");
    enm_non_trimmed_quad_rule_settings
        .value("integration_method", NonTrimmedQuadratureRuleSettings::integration_method)
    ;
    py::implicitly_convertible<std::string, NonTrimmedQuadratureRuleSettings>();

    py::enum_<ConditionSettings> enm_condition_settings(m, "ConditionSettings");
    enm_condition_settings
        .value("condition_id", ConditionSettings::condition_id)
        .value("condition_type", ConditionSettings::condition_type)
        .value("input_filename", ConditionSettings::input_filename)
        .value("modulus", ConditionSettings::modulus)
        .value("value", ConditionSettings::value)
        .value("penalty_factor", ConditionSettings::penalty_factor)
    ;
    py::implicitly_convertible<std::string, ConditionSettings>();

    /// Export Dictionary
    typedef py::class_<SettingsBaseType> SettingsBaseTypeBinderType;
    py::class_<SettingsBaseType> settings_base_type_binder(m,"Dictionary");

    DictionaryBinderHelper<SettingsBaseTypeBinderType, Main>(settings_base_type_binder);
    DictionaryBinderHelper<SettingsBaseTypeBinderType, GeneralSettings>(settings_base_type_binder);
    DictionaryBinderHelper<SettingsBaseTypeBinderType, MeshSettings>(settings_base_type_binder);
    DictionaryBinderHelper<SettingsBaseTypeBinderType, TrimmedQuadratureRuleSettings>(settings_base_type_binder);
    DictionaryBinderHelper<SettingsBaseTypeBinderType, NonTrimmedQuadratureRuleSettings>(settings_base_type_binder);
    DictionaryBinderHelper<SettingsBaseTypeBinderType, ConditionSettings>(settings_base_type_binder);

    /// Export Parameters
    py::class_<Settings, SettingsBaseType>(m,"Settings")
        .def(py::init<>())
        .def("CreateNewConditionSettings", &Settings::CreateNewConditionSettings, py::return_value_policy::reference_internal)
    ;

    }

} // End namespace Python
} // End namespace queso

