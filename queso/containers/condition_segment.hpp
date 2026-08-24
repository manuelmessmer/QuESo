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

#pragma once

//// STL includes
#include <optional>

//// Project includes
#include "queso/containers/triangle_mesh.hpp"
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

/// @class  ConditionSegment
/// @author Manuel Messmer
/// @brief  Segment of a condition clipped to the background-grid element boundaries.
///         Stores the parent element id and, if the parent element is active, a reference to it.
///         Additionally stores the clipped section of the triangle mesh.
/// @todo   Add boundary integration points that can eventually be used by a boundary moment-fitting scheme.
template<typename TElementViewType>
class ConditionSegment
{
public:
    ///@name Type definitions
    ///@{
    using ElementViewType = TElementViewType;
    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor
    /// @param Index index in the background grid. @see GridIndexer.
    /// @param rTriangleMesh Triangle mesh of this segment. ConditionSegment takes ownership.
    ConditionSegment(IndexType Index, TriangleMesh&& rTriangleMesh)
        : mParentCellIndex(Index), mTriangleMesh(std::move(rTriangleMesh))
    {}

    /// @brief Constructor
    /// @param Index index in the background grid. @see GridIndexer.
    /// @param rTriangleMesh Triangle mesh of this segment. ConditionSegment takes ownership.
    /// @param rElement parent element.
    ConditionSegment(IndexType Index, TriangleMesh&& rTriangleMesh, const ElementViewType& rElement)
        : mParentCellIndex(Index), mTriangleMesh(std::move(rTriangleMesh)), mParentElement(rElement)
    {}

    /// Destructor
    ~ConditionSegment() = default;
    /// Copy constructor
    ConditionSegment(const ConditionSegment& rOther) = delete;
    /// Assignment operator
    ConditionSegment& operator=(const ConditionSegment& rOther) = delete;
    /// Move constructor
    ConditionSegment(ConditionSegment&& rOther) noexcept = default;
    /// Move assignment operator
    ConditionSegment& operator=(ConditionSegment&& rOther) noexcept = default;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns the triangle mesh representing the condition segment.
    /// @return const TriangleMesh&
    [[nodiscard]] const TriangleMesh& GetTriangleMesh() const
    { return mTriangleMesh; }

    /// @brief Returns the parent background-grid cell index.
    /// @return Index of the parent cell.
    [[nodiscard]] IndexType GetParentCellIndex() const noexcept
    { return mParentCellIndex; }

    /// @brief Returns true if the condition segment is contained within an active parent element.
    /// @return bool
    [[nodiscard]] bool IsInActiveElement() const noexcept
    { return mParentElement.has_value(); }

    /// @brief Returns the optional parent element.
    /// @return const optional parent-element reference.
    [[nodiscard]] const std::optional<ElementViewType>& GetParentElement() const noexcept
    { return mParentElement; }

private:
    ///@}
    ///@name Private member variables
    ///@{

    IndexType mParentCellIndex{};
    TriangleMesh mTriangleMesh{};
    std::optional<ElementViewType> mParentElement{};
    // BoundaryIntegrationPointVectorType mIntegrationPoints;

    ///@}
};  // End class ConditionSegment
///@}
}  // namespace queso
