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

#ifndef FLOOD_FILL_TESTER_INCLUDE_HPP
#define FLOOD_FILL_TESTER_INCLUDE_HPP

#include "queso/embedding/flood_fill.h"

namespace queso {
namespace Testing {

// Make protected funtions public for testing.
// This class shall only be used for testing!
class FloodFillTester : public FloodFill {
public:
    FloodFillTester(BRepOperator* pBrepOperator, const Parameters& rParameters)
        : FloodFill(pBrepOperator, rParameters)
    {
    }

    using FloodFill::ClassifyElementsForTest;
};

} // End namespace Tester
} // End namespace queso

#endif // FLOOD_FILL_TESTER_INCLUDE_HPP