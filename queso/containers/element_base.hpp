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

//// Project includes
#include "queso/containers/cell_domain.hpp"
#include "queso/containers/element_core.hpp"
#include "queso/containers/element_view.hpp"
#include "queso/utilities/mapping_utilities.hpp"

namespace queso {

namespace concepts {
    template<typename T>
    concept CellLikeDomain = std::same_as<T, CellDomain>;
}
///@name QuESo Classes
///@{

/// @brief Shared implementation for trimmed and untrimmed element types.
/// @details Stores the common element core, exposes integration-point access in
/// global or parametric space, and forwards active-domain queries to the stored
/// domain object.
/// @tparam TIntegrationPoint Volume integration-point type.
/// @tparam TBoundaryIntegrationPoint Boundary integration-point type.
/// @tparam TDomain Active-domain type.
/// @tparam TIsTrimmed Compile-time trimmed flag exposed by IsTrimmed().
template<
    concepts::IntegrationPoint TIntegrationPoint,
    concepts::BoundaryIntegrationPoint TBoundaryIntegrationPoint,
    typename TDomain,
    bool TIsTrimmed>
class ElementBase
{
public:
    ///@name Type definitions
    ///@{

    using CoreType = ElementCore<TIntegrationPoint, TBoundaryIntegrationPoint>;
    using ElementViewType = ElementView<TIntegrationPoint, TBoundaryIntegrationPoint>;
    using IntegrationPointType = typename CoreType::IntegrationPointType;
    using BoundaryIntegrationPointType = typename CoreType::BoundaryIntegrationPointType;
    using IntegrationPointVectorType = typename CoreType::IntegrationPointVectorType;
    using BoundaryIntegrationPointVectorType = std::vector<BoundaryIntegrationPointType>;
    using DomainType = TDomain;
    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructs the shared element base from its core data and domain.
    /// @param rCore Element core with id, bounds, mapper, and integration points.
    /// @param rDomain Active-domain representation owned by the element.
    ElementBase(CoreType&& rCore, DomainType&& rDomain) : mCore(std::move(rCore)), mDomain(std::move(rDomain))
    {}

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns a lightweight read-only view of this element.
    /// @return ElementViewType
    [[nodiscard]] ElementViewType View() const noexcept
    { return ElementViewType(mCore, TIsTrimmed); }

    /// @brief Returns the element id.
    /// @return IndexType
    [[nodiscard]] IndexType GetId() const noexcept
    { return mCore.id; }

    /// @brief Returns the element bounds in the requested coordinate space.
    /// @tparam TSpace Coordinate space.
    /// @return const BoundingBoxType&
    template<CoordinateSpace TSpace>
    [[nodiscard]] const BoundingBoxType& GetCellBounds() const noexcept
    {
        if constexpr (TSpace == CoordinateSpace::global) {
            return mCore.bounds.global;
        } else {
            return mCore.bounds.parametric;
        }
    }

    /// @brief Returns mutable access to the stored parametric integration points.
    /// @return IntegrationPointVectorType&
    [[nodiscard]] IntegrationPointVectorType& GetIntegrationPoints() noexcept
    { return mCore.integration_points; }

    /// @brief Returns integration points in the requested coordinate space.
    /// @details Parametric access returns the stored integration-point container.
    /// Global access returns a lazily transformed view using the element-local
    /// mapping and Jacobian determinant.
    /// @tparam TSpace Coordinate space.
    /// @return Container reference or transformed view depending on TSpace.
    template<CoordinateSpace TSpace = CoordinateSpace::parametric>
    [[nodiscard]] decltype(auto) GetIntegrationPoints() const
    {
        if constexpr (TSpace == CoordinateSpace::parametric) {
            return static_cast<const IntegrationPointVectorType&>(mCore.integration_points);
        } else {
            static_assert(
                std::constructible_from<IntegrationPointType, PointType, double>,
                "Global integration-point access requires IntegrationPointType(PointType, double)."
            );

            const double det_j = DetJ();
            return mCore.integration_points | std::views::transform([this, det_j](const auto& rPoint) {
                       return IntegrationPointType(
                           mapping::ToGlobal(rPoint.Point(), mCore.bounds), rPoint.Weight() * det_j
                       );
                   });
        }
    }

    /// @brief Returns the determinant of the element-local mapping Jacobian.
    /// @return double
    [[nodiscard]] double DetJ() const noexcept
    { return mapping::DetJ(mCore.bounds); }

    /// @brief Returns whether this element represents a trimmed domain.
    /// @return bool
    [[nodiscard]] bool IsTrimmed() const noexcept
    { return TIsTrimmed; }

    /// @brief Returns the active-domain bounds in the requested space.
    /// @tparam TSpace Coordinate space.
    /// @return BoundingBoxType
    template<CoordinateSpace TSpace>
    [[nodiscard]] BoundingBoxType GetActiveDomainBounds() const
    { return mDomain.template GetBounds<TSpace>(mCore.bounds); }

    /// @brief Returns the boundary mesh of the active domain.
    /// @return TriangleMeshView
    [[nodiscard]] TriangleMeshView GetActiveDomainBoundaryMesh() const
    {
        if constexpr (concepts::CellLikeDomain<DomainType>) {
            return mDomain.GetBoundaryMesh(mCore.bounds);
        } else {
            return mDomain.GetBoundaryMesh();
        }
    }

    /// @brief Returns boundary integration points of the active domain.
    /// @tparam TBoundaryIpType Boundary integration-point type.
    /// @tparam TSpace Coordinate space of the returned points.
    /// @return std::vector<TBoundaryIpType>
    template<typename TBoundaryIpType = BoundaryIntegrationPointType, CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] std::vector<TBoundaryIpType> GetActiveDomainBoundaryIps() const
    { return mDomain.template GetBoundaryIps<TBoundaryIpType, TSpace>(mCore.bounds); }

    /// @brief Returns whether a point lies inside the active domain.
    /// @tparam TSpace Coordinate space of the query point.
    /// @param rPoint Query point.
    /// @return bool
    template<CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] bool IsInsideActiveDomain(const PointType& rPoint) const
    { return mDomain.template IsInside<TSpace>(rPoint, mCore.bounds); }

    /// @brief Returns the intersection state of a query box against the active domain.
    /// @tparam TSpace Coordinate space of the query bounds.
    /// @param rLowerBound Lower bound of the query box.
    /// @param rUpperBound Upper bound of the query box.
    /// @param Tolerance Shrink tolerance applied to the query box.
    /// @return IntersectionStateType
    template<CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] IntersectionStateType GetIntersectionState(
        const PointType& rLowerBound,
        const PointType& rUpperBound,
        double Tolerance = SNAPTOL
    ) const
    { return mDomain.template GetIntersectionState<TSpace>(rLowerBound, rUpperBound, mCore.bounds, Tolerance); }

    ///@}
protected:
    ///@name Protected members
    ///@{

    CoreType mCore;
    DomainType mDomain;
    ///@}
};
///@}
}// namespace queso
