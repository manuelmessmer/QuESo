// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

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