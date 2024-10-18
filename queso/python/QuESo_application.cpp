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

//// Project includes
#include "queso/python/define_python.hpp"
// To export
#include "queso/python/add_settings_to_python.h"
#include "queso/python/add_model_info_to_python.h"
#include "queso/python/add_containers_to_python.h"
#include "queso/python/add_globals_to_python.h"
#include "queso/python/add_utilities_to_python.h"
#include "queso/QuESo_main.h"

namespace queso {
namespace Python {
namespace py = pybind11;

PYBIND11_MODULE(QuESo_Application,m) {

    m.doc() = "This is a Python binding for QuESo";

    m.def("PrintLogo", []()
    {
    QuESo_INFO << " Importing QuESo \n"
        << "   ____        ______  _____        \n"
        << "  / __ \\      |  ____|/ ____|       \n"
        << " | |  | |_   _| |__  | (___   ___   \n"
        << " | |  | | | | |  __|  \\___ \\ / _ \\  \n"
        << " | |__| | |_| | |____ ____) | (_) | \n"
        << "  \\___\\_\\\\__,_|______|_____/ \\___/  \n"
        << "\t Quadrature for Embedded Solids \n\n";
    }, "Print Logo");

    AddGlobalsToPython(m);
    AddSettingsToPython(m);
    AddModelInfoToPython(m);
    AddContainersToPython(m);
    AddUtilitiesToPython(m);

}

}// End namespace Python
}// End namespace queso
