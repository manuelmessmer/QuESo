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
#include "queso/python/add_model_info_to_python.h"
#include "queso/python/dictionary_binder_helper.hpp"
// To export
#include "queso/includes/model_info.hpp"

// Note: PYBIND11_MAKE_OPAQUE can not be captured within namespace
typedef std::vector<queso::ModelInfoBaseType> ModelInfoDictionaryVectorType;
PYBIND11_MAKE_OPAQUE(ModelInfoDictionaryVectorType);

namespace queso {
namespace Python {

namespace py = pybind11;

void AddModelInfoToPython(pybind11::module& m) {

    /// Export ModelInfoDictionary
    py::class_<ModelInfoBaseType> model_info_base_type_binder(m,"ModelInfoDictionary");

    DictionaryBinderHelper<ModelInfoBaseType>(model_info_base_type_binder);

    /// Export ModelInfoDictionaryVectorType (Required for list of subdicionaries)
    py::bind_vector<ModelInfoDictionaryVectorType, Unique<ModelInfoDictionaryVectorType>>
        (m, "ModelInfoDictionaryVector")
    ;

    /// Export ModelInfo
    // Note: ModelInfo is always created within C++. Therefore, no constructor required.
    py::class_<ModelInfo, ModelInfoBaseType>(m,"ModelInfo")
    ;

} // End AddModelInfoToPython

} // End namespace Python
} // End namespace queso

