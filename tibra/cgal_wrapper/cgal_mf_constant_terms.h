// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CGAL_MF_CONSTANT_TERMS_H
#define CGAL_MF_CONSTANT_TERMS_H

// External includes
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <array>

// Project includes
#include "utilities/parameters.h"
#include "geometries/element.h"

namespace cgal {

///@name TIBRA Classes
///@{

/**
 * @class  CGAL constant terms
 * @author Manuel Messmer
 * @brief Provides functions to compute the constant terms of the moment fitting equation.
*/
class ConstantTerms{

public:
    ///@name Type Definitions
    ///@{
    typedef Element::IntegrationPointVectorType IntegrationPointVectorType;
    typedef boost::numeric::ublas::vector<double> VectorType;

    ///@}
    ///@name Operatios
    ///@{

    ///@brief Construct constant terms for moment fitting equation.
    ///@param rElement
    ///@param [out] rConstantTerms
    ///@param rParam
    static void Compute(const Element& rElement, VectorType& rConstantTerms, const Parameters& rParam);

    ///@}

}; // End Class
///@}

} // End cgal namespace

#endif // CGAL_MF_CONSTANT_TERMS_H