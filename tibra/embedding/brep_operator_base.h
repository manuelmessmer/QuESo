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
    typedef std::vector<IntersectionStatusType> StatusVectorType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    BRepOperatorBase(const Parameters& rParameters) : mParameters(rParameters)
    {
        TIBRA_ERROR_IF("BRepOperatorBase::Constructor", !rParameters.Has<PointType>("lower_bound") ) << "Parameter does not contain 'lower_bound'.\n";
        TIBRA_ERROR_IF("BRepOperatorBase::Constructor", !rParameters.Has<PointType>("upper_bound") ) << "Parameter does not contain 'upper_bound'.\n";
        TIBRA_ERROR_IF("BRepOperatorBase::Constructor", !rParameters.Has<Vector3i>("number_of_elements") ) << "Parameter does not contain 'number_of_elements'.\n";
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
    ///@note Calls: GetIntersectionState(const PointType& rLowerBound,  const PointType& rUpperBound, double Tolerance = SNAPTOL)
    ///@param rElement
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    virtual IntersectionStatus GetIntersectionState(const Element& rElement){
        const auto& lower_bound = rElement.GetLowerBound();
        const auto& upper_bound = rElement.GetUpperBound();
        return GetIntersectionState(lower_bound, upper_bound);
    }

    ///@brief Returns intersections state of AABB.
    ///@param rLowerBound Lower bound of AABB.
    ///@param rUpperBound Lower bound of AABB.
    ///@param Tolerance Default is SNAPTOL. Slightly reduces aabb such that "touches" are not considered as trimming.
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    virtual IntersectionStatus GetIntersectionState(const PointType& rLowerBound,  const PointType& rUpperBound, double Tolerance = SNAPTOL) const = 0;

    /// @brief Returns a ptr to a vector that holds the states of each element. Vector is ordered according to index -> see: Mapper.
    /// @brief This function runs a flood fill repeatively and classifies each group based on the bounding elements that are trimmed. Each element that borders a trimmed
    ///        element is tested via local ray tracing and marked as inside or outside. The majority vote decides about the classification of each group.
    /// @return Unique<StatusVectorType>.
    virtual Unique<StatusVectorType> GetElementClassifications() const = 0;

    /// @brief Returns true, if AABB is intersected by at least one triangle.
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @param Tolerance Reduces size of AABB.
    /// @return bool.
    virtual bool IsTrimmed(const PointType& rLowerBound,  const PointType& rUpperBound, double Tolerance = SNAPTOL) const = 0;

    /// @brief Returns ptr to trimmed domain.
    /// @param rLowerBound Lower bound of AABB.
    /// @param rUpperBound Lower bound of AABB.
    /// @param rParam Parameterss
    /// @return TrimmedDomainBasePtrType (Unique)
    virtual TrimmedDomainBasePtrType pGetTrimmedDomain(const PointType& rLowerBound, const PointType& rUpperBound ) const = 0;

    /// @brief Returns true if rPoint lies on bounded side of clipped mesh (clipped by AABB).
    ///        Ray tracing through the center of at least 10 triangles (or maximum number of triangles, if n_max < 10) is performed.
    ///        The majority decides about the classification of rPoint. Note that this function is much more efficient than IsInside.
    ///        However, rPoint must be close to AABB. This is e.g. used to classify an aabb next to a trimmed aabb (see: FloodFlow()).
    /// @param rPoint Query Point.
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @return bool
    virtual bool OnBoundedSideOfClippedSection( const PointType& rPoint, const PointType& rLowerBound, const PointType& rUpperBound ) const = 0;

protected:
    ///@}
    ///@name Protected memmbers
    ///@{

    const Parameters& mParameters;

    ///@}
}; // End BRepOperatorBase class
///@} // End TIBRA classes

} // End namespace tibra

#endif // BREP_OPERATOR_BASE_INCLUDE_H