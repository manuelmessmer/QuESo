// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MOMENT_FITTING_UTILITIES_INCLUDE_H
#define MOMENT_FITTING_UTILITIES_INCLUDE_H

// External includes
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <array>
#include <variant>

// Project includes
#include "geometries/element.h"
#include "utilities/parameters.h"



class MomentFitting{
public:
    typedef Element::IntegrationPointVectorType IntegrationPointVectorType;
    typedef boost::numeric::ublas::vector<double> VectorType;

    static void CreateIntegrationPointsTrimmed(Element& rElement, const Parameters& rParam);

private:
    static double CreateIntegrationPointsTrimmed(Element& rElement, const VectorType& rConstantTerms, const int PointDistributionFactor, const Parameters& rParam);

    static void DistributeInitialIntegrationPoints(const Element& rElement, IntegrationPointVectorType& rIntegrationPoint, const int PointDistributionFactor, const Parameters& rParam);

    static void ComputeConstantTerms(const Element& rElement, VectorType& rConstantTerms, const Parameters& rParam);
}; // End Class

#endif // MOMENT_FITTING_UTILITIES_INCLUDE_H