// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIMMED_ELEMENT_TESTER_INCLUDE_HPP
#define TRIMMED_ELEMENT_TESTER_INCLUDE_HPP

namespace tibra {
namespace Testing {

// Make protected funtions public for testing.
// This function shall only be used for testing!
class QuadratureTrimmedElementTester : public QuadratureTrimmedElement {
public:
    using QuadratureTrimmedElement::DistributeIntegrationPoints;
    using QuadratureTrimmedElement::ComputeConstantTerms;
    using QuadratureTrimmedElement::MomentFitting;
    using QuadratureTrimmedElement::PointElimination;
};

} // End namespace Tester
} // End namespace tibra

#endif // TRIMMED_ELEMENT_TESTER_INCLUDE_HPP