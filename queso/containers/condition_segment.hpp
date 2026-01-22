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
#include <vector>

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/containers/triangle_mesh_interface.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  ConditionSegment
 * @author Manuel Messmer
 * @brief  Segment of condition that is clipped to the element boundaries of the background grid.
 *         Stores the index of the parent element (see: GridIndexer). If parent element is active, also stores a ptr to the parent element.
 *         Additionally, stores the clipped section of the triangle mesh, and the corresponding boundary integration points (not yet).
 * @todo   Add BoundaryIntegrationPoints that can eventually be used by a boundary MomentFitting scheme.
**/
template<typename TElementType>
class ConditionSegment {
public:

    ///@name Type Definitions
    ///@{
    typedef TElementType ElementType;
    typedef typename ElementType::BoundaryIntegrationPointType BoundaryIntegrationPointType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param Index of BackgroundGrid @see GridIndexer.
    /// @param pTriangleMesh ptr to triangle mesh of this segment.
    ConditionSegment(IndexType Index, Unique<TriangleMeshInterface>& pTriangleMesh)
        : mBackgroundGridIndex(Index), mpParentElement(nullptr), mpTriangleMesh(std::move(pTriangleMesh))
    {
    }

    /// @brief Constructor
    /// @param Index of BackgroundGrid @see GridIndexer.
    /// @param pElement Ptr to parent element.
    /// @param pTriangleMesh ptr to triangle mesh of this segment.
    ConditionSegment(IndexType Index, const ElementType* pElement, Unique<TriangleMeshInterface>& pTriangleMesh)
        : mBackgroundGridIndex(Index), mpParentElement(pElement), mpTriangleMesh(std::move(pTriangleMesh))
    {
    }

    // Destructor
    ~ConditionSegment() = default;
    // Copy constructor
    ConditionSegment(const ConditionSegment& rOther) = delete;
    // Assignement operator
    ConditionSegment& operator=(const ConditionSegment& rOther) = delete;
    /// Move constructor
    ConditionSegment(ConditionSegment&& rOther) noexcept = default;
    /// Move assignement operator
    ConditionSegment& operator=(ConditionSegment&& rOther) noexcept = delete;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Return const ref to the triangle mesh, representing the condition segment.
    /// @return  const TriangleMeshInterface&
    const TriangleMeshInterface& GetTriangleMesh() const {
        return *mpTriangleMesh;
    }

    /// @brief Returns true if ConditionSegment is contained within in active parent element.
    /// @return bool
    bool IsInActiveElement() const {
        return (mpParentElement) != 0;
    }

private:

    ///@}
    ///@name Private member variables
    ///@{

    const IndexType mBackgroundGridIndex;
    const ElementType* mpParentElement;
    Unique<const TriangleMeshInterface> mpTriangleMesh;

    ///@}
}; // End class ConditionSegment
///@}
} // End queso namespace.

#endif // End CONDITION_SEGMENT_INCLUDE_HPP


