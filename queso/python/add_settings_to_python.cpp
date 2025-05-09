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
#include "queso/python/dictionary_binder_helper.hpp"
// To export
#include "queso/includes/settings.hpp"


// Note: PYBIND11_MAKE_OPAQUE can not be captured within namespace
using SettingsDictionaryVectorType = std::vector<queso::SettingsBaseType>;
PYBIND11_MAKE_OPAQUE(SettingsDictionaryVectorType);

namespace queso {
namespace Python {

namespace py = pybind11;

void AddSettingsToPython(pybind11::module& m) {

    /// Export Dictionary
    py::class_<SettingsBaseType> settings_base_type_binder(m,"SettingsDictionary");

    DictionaryBinderHelper<SettingsBaseType>(settings_base_type_binder);

    /// Export SettingsDictionaryVectorType (Required for list of subdicionaries)
    py::bind_vector<SettingsDictionaryVectorType>(m, "SettingsDictionaryVector");

    /// Export Settings
    py::class_<Settings, SettingsBaseType>(m,"Settings")
        .def(py::init<>())
        .def("CreateNewConditionSettings", &Settings::CreateNewConditionSettings, py::return_value_policy::reference_internal)
    ;

} // End AddSettingsToPython

} // End namespace Python
} // End namespace queso

