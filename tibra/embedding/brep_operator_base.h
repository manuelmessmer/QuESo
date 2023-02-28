// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BREP_OPERATOR_BASE_INCLUDE_H
#define BREP_OPERATOR_BASE_INCLUDE_H

//// STL includes
#include <memory>
//// Project includes
#include "containers/element.hpp"
#include "embedding/trimmed_domain_base.h"
#include "utilities/parameters.h"
#include "define.hpp"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  BRepOperatorBase
 * @author Manuel Messmer
 * @brief  Base class for BRepOperator
*/
class BRepOperatorBase {

public:
    ///@name Type Definitions
    ///@{
    typedef Unique<TrimmedDomainBase> TrimmedDomainBasePtrType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    BRepOperatorBase(const Parameters& rParameters) : mParameters(rParameters)
    {
    }

    /// Destructor
    virtual ~BRepOperatorBase() {}

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TriangleMesh.
    ///@param rPoint
    ///@return bool
    virtual bool IsInside(const PointType& rPoint) const = 0;

    ///@brief Returns intersections state of element.
    ///@param rElement
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    virtual IntersectionStatus GetIntersectionState(const Element& rElement){
        const auto& lower_bound = rElement.GetLowerBound();
        const auto& upper_bound = rElement.GetUpperBound();
        return GetIntersectionState(lower_bound, upper_bound);
    }

    ///@brief Returns intersections state of AABB.
    ///@param rLowerBound of AABB.
    ///@param rUpperBound of AABB.
    ///@param Tolerance Default is SNAPTOL. Slightly reduces aabb such that "touches" are not considered as trimming.
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    virtual IntersectionStatus GetIntersectionState(const PointType& rLowerBound,  const PointType& rUpperBound, double Tolerance = SNAPTOL) const = 0;

    /// @brief Returns ptr to trimmed domain.
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @param rParam Parameterss
    /// @return TrimmedDomainBasePtrType (Unique)
    virtual TrimmedDomainBasePtrType pGetTrimmedDomain(const PointType& rLowerBound, const PointType& rUpperBound ) const = 0;

protected:
    const Parameters& mParameters;
    ///@}
}; // End BRepOperatorBase class
///@} // End TIBRA classes

} // End namespace tibra

#endif // BREP_OPERATOR_BASE_INCLUDE_H