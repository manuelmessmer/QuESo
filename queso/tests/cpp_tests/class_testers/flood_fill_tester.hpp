// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef FLOOD_FILL_TESTER_INCLUDE_HPP
#define FLOOD_FILL_TESTER_INCLUDE_HPP

#include "embedding/flood_fill.h"

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