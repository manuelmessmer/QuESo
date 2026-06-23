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
#include "queso/containers/element_base.hpp"
#include "queso/embedding/trimmed_domain.h"

namespace queso {

///@name QuESo Classes
///@{

/// @class TrimmedElement
/// @author Manuel Messmer
/// @brief Element representation for trimmed active domains.
/// @details Extends ElementBase with a TrimmedDomain and exposes the common
/// element API for embedded geometries. The element stores integration points
/// in parametric space and delegates active-domain queries to the underlying
/// trimmed-domain representation.
/// @tparam TIntegrationPointType Volume integration-point type.
/// @tparam TBoundaryIntegrationPointType Boundary integration-point type.
template<
    concepts::IntegrationPoint TIntegrationPointType,
    concepts::BoundaryIntegrationPoint TBoundaryIntegrationPointType>
class TrimmedElement : public ElementBase<TIntegrationPointType, TBoundaryIntegrationPointType, TrimmedDomain, true>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ElementBase<TIntegrationPointType, TBoundaryIntegrationPointType, TrimmedDomain, true>;

    using CoreType = typename BaseType::CoreType;
    using ElementViewType = typename BaseType::ElementViewType;
    using IntegrationPointType = typename BaseType::IntegrationPointType;
    using BoundaryIntegrationPointType = typename BaseType::BoundaryIntegrationPointType;
    using IntegrationPointVectorType = typename BaseType::IntegrationPointVectorType;
    using BoundaryIntegrationPointVectorType = typename BaseType::BoundaryIntegrationPointVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructs a trimmed element from its components.
    /// @param Id Unique element identifier.
    /// @param rBounds Element bounds in global and parametric space.
    /// @param rDomain Trimmed-domain representation.
    /// @param IntegrationPoints Volume integration points in parametric space.
    TrimmedElement(
        IndexType Id,
        const ElementBounds& rBounds,
        TrimmedDomain&& rDomain,
        IntegrationPointVectorType IntegrationPoints = {}
    )
        : TrimmedElement(CoreType{ Id, rBounds, std::move(IntegrationPoints) }, std::move(rDomain))
    {}

    /// @brief Constructs a trimmed element from a preassembled core and domain.
    /// @param rCore Element core data.
    /// @param rDomain Trimmed-domain representation.
    TrimmedElement(CoreType&& rCore, TrimmedDomain&& rDomain) : BaseType(std::move(rCore), std::move(rDomain))
    {}

    ///@}
};

///@}
}// namespace queso
