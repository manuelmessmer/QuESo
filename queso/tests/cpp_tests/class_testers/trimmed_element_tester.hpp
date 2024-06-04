//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef TRIMMED_ELEMENT_TESTER_INCLUDE_HPP
#define TRIMMED_ELEMENT_TESTER_INCLUDE_HPP

#include "queso/quadrature/trimmed_element.hpp"

namespace queso {
namespace Testing {

// Make protected funtions public for testing.
// This function shall only be used for testing!
template<typename TElementType>
class QuadratureTrimmedElementTester : public QuadratureTrimmedElement<TElementType> {
public:
    using QuadratureTrimmedElement<TElementType>::DistributeIntegrationPoints;
    using QuadratureTrimmedElement<TElementType>::ComputeConstantTerms;
    using QuadratureTrimmedElement<TElementType>::MomentFitting;
    using QuadratureTrimmedElement<TElementType>::PointElimination;
};

} // End namespace Tester
} // End namespace queso

#endif // TRIMMED_ELEMENT_TESTER_INCLUDE_HPP