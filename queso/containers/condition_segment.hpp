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

#ifndef CONDITION_SEGMENT_INCLUDE_HPP
#define CONDITION_SEGMENT_INCLUDE_HPP

//// STL includes
#include "queso/includes/settings.hpp"
#include "queso/containers/triangle_mesh.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  ConditionSegment
 * @author Manuel Messmer
 * @brief  Segment of condition that is clipped to the element boundaries of the background grid.
 *         Stores a ptr to the parent element, the clipped section of the triangle mesh, and the corresponding boundary integration points.
**/
template<typename TElementType>
class ConditionSegment {
public:

    ///@name Type Definitions
    ///@{
    typedef TElementType ElementType;
    typedef typename ElementType::BoundaryIntegrationPointType BoundaryIntegrationPointType;
    typedef std::vector<BoundaryIntegrationPointType> BoundaryIntegrationPointVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param pElement Ptr to parent element.
    /// @param pTriangleMesh ptr to triangle mesh of this segment.
    ConditionSegment(const ElementType* pElement, Unique<TriangleMeshInterface>& pTriangleMesh)
        : mpParentElement(pElement), mpTriangleMesh(std::move(pTriangleMesh))
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Return const ref to the triangle mesh, representing the boundary.
    /// @return  const TriangleMeshInterface&
    const TriangleMeshInterface& GetTriangleMesh() const {
        return *mpTriangleMesh;
    }

    /// @brief Returns true if ConditionSegment is contained within in active parent element.
    /// @return bool
    bool IsInActiveElement() const {
        return (mpParentElement) ? true : false;
    }

    // Unique<std::vector<BoundaryIntegrationPointType>> pGetIntegrationPoints() const {
    //     // Pointer to boundary integration points
    //     auto p_boundary_ips = MakeUnique<BoundaryIntegrationPointVectorType>();

    //     p_boundary_ips->reserve(mpTriangleMesh->NumOfTriangles()*12UL);
    //     for( IndexType triangle_id = 0; triangle_id < mpTriangleMesh->NumOfTriangles(); ++triangle_id ){
    //         IndexType method = 0;
    //         auto p_new_points = mpTriangleMesh->pGetIPsGlobal<BoundaryIntegrationPointType>(triangle_id, method);
    //         p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
    //     }

    //     return p_boundary_ips;
    // }

private:

    ///@}
    ///@name Private member variables
    ///@{

    const ElementType* mpParentElement;
    const Unique<TriangleMeshInterface> mpTriangleMesh;
    BoundaryIntegrationPointVectorType mIntegrationPoints;

    ///@}
}; // End class ConditionSegment
///@}
} // End queso namespace.

#endif // End CONDITION_SEGMENT_INCLUDE_HPP


