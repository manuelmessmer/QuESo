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
#include <functional>

//// Project includes
#include "queso/containers/element_core.hpp"
#include "queso/utilities/mapping_utilities.hpp"

namespace queso {

///@name QuESo classes
///@{

/// @class  ElementView
/// @author Manuel Messmer
/// @brief  Lightweight, non-owning read-only view of an element.
/// @details Holds a reference_wrapper to the immutable ElementCore and caches the
/// element-local Jacobian determinant computed at construction time.  The view is
/// cheap to copy and suitable for passing across API boundaries without transferring
/// ownership.
/// @tparam TIntegrationPoint Volume integration-point type.
/// @tparam TBoundaryIntegrationPoint Boundary integration-point type.
template<concepts::IntegrationPoint TIntegrationPoint, concepts::BoundaryIntegrationPoint TBoundaryIntegrationPoint>
class ElementView
{
public:
    ///@name Type Definitions
    ///@{

    using CoreType = ElementCore<TIntegrationPoint, TBoundaryIntegrationPoint>;
    using IntegrationPointType = typename CoreType::IntegrationPointType;
    using BoundaryIntegrationPointType = typename CoreType::BoundaryIntegrationPointType;
    using IntegrationPointVectorType = typename CoreType::IntegrationPointVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructs a view from an element core.
    /// @details The Jacobian determinant is computed once from the element-local
    /// bounds stored in the core and cached for the lifetime of the view.
    /// @param rCore  Element core to observe.  Must outlive this view.
    /// @param IsTrimmed  Whether the element represents a trimmed domain.
    ElementView(const CoreType& rCore, bool IsTrimmed) noexcept
        : mCore(rCore), mDetJ(mapping::DetJ(rCore.bounds)), mIsTrimmed(IsTrimmed)
    {}

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns the element id.
    /// @return IndexType
    [[nodiscard]] IndexType GetId() const noexcept
    { return mCore.get().id; }

    /// @brief Returns the element bounds in the requested coordinate space.
    /// @tparam TSpace Coordinate space.
    /// @return const BoundingBoxType&
    template<CoordinateSpace TSpace>
    [[nodiscard]] const BoundingBoxType& GetCellBounds() const noexcept
    {
        if constexpr (TSpace == CoordinateSpace::global) {
            return mCore.get().bounds.global;
        } else {
            return mCore.get().bounds.parametric;
        }
    }

    /// @brief Returns integration points in the requested coordinate space.
    /// @details Parametric access returns the stored integration-point container
    /// directly.  Global access returns a lazily transformed view where each point is
    /// mapped to global coordinates and its weight is scaled by DetJ().
    /// @tparam TSpace Coordinate space.
    /// @return Container reference (parametric) or transformed view (global).
    template<CoordinateSpace TSpace = CoordinateSpace::parametric>
    [[nodiscard]] decltype(auto) GetIntegrationPoints() const
    {
        if constexpr (TSpace == CoordinateSpace::parametric) {
            return static_cast<const IntegrationPointVectorType&>(mCore.get().integration_points);
        } else {
            static_assert(
                std::constructible_from<IntegrationPointType, PointType, double>,
                "Global integration-point access requires IntegrationPointType(PointType, double)."
            );
            const auto& r_core = mCore.get();
            return r_core.integration_points | std::views::transform([this, &r_core](const auto& rPoint) {
                       return IntegrationPointType(
                           mapping::ToGlobal(rPoint.Point(), r_core.bounds), rPoint.Weight() * mDetJ
                       );
                   });
        }
    }

    /// @brief Returns the cached element-local Jacobian determinant.
    /// @return double
    [[nodiscard]] double DetJ() const noexcept
    { return mDetJ; }

    /// @brief Returns whether this view represents a trimmed element.
    /// @return bool
    [[nodiscard]] bool IsTrimmed() const noexcept
    { return mIsTrimmed; }

    ///@}
private:
    ///@name Private Members
    ///@{

    std::reference_wrapper<const CoreType> mCore;
    double mDetJ{};
    bool mIsTrimmed{};

    ///@}
};
///@}
}// namespace queso
