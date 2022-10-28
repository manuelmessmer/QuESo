// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MF_CONSTANT_TERMS_H
#define MF_CONSTANT_TERMS_H

// External includes
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <array>

// Project includes
#include "utilities/parameters.h"
#include "geometries/element.h"
#include "geometries/boundary_integration_point.h"


///@name TIBRA Classes
///@{

/**
 * @class  Constant terms
 * @author Manuel Messmer
 * @brief New version of costant terms. Will be used once all CGAl dependencies are removed.
*/
class ConstantTerms{

public:
    ///@name Type Definitions
    ///@{
    typedef std::vector<BoundaryIntegrationPoint> BoundaryIPsVectorType;
    typedef std::unique_ptr<BoundaryIPsVectorType> BoundaryIPsVectorPtrType;
    typedef boost::numeric::ublas::vector<double> VectorType;
    ///@}

    ///@name Operations
    ///@{

    ///@brief Compute constant terms of moment fitting equation.
    ///@param rBoundaryIps
    ///@param rElement
    ///@param rConstantTerms
    ///@param rParam
    ///@todo Remove element and pass lower_bound and upper_bound.
    static void Compute(BoundaryIPsVectorPtrType& rBoundaryIps, const Element& rElement, VectorType& rConstantTerms, const Parameters& rParam);

    ///@}

}; // End Class

///@}

#endif // MF_CONSTANT_TERMS_H