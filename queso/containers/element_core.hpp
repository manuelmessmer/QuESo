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
#include <vector>

//// Project includes
#include "queso/includes/define.hpp"

namespace queso {

/// @brief Affine element bounds in global and parametric space.
struct ElementBounds
{
    BoundingBoxType global;
    BoundingBoxType parametric;
};

///@name QuESo Classes
///@{

/// @class ElementCore
/// @author Manuel Messmer
/// @brief Shared immutable data stored by all element types.
/// @details Holds the element identifier, element-local bounds, and volume
/// integration points. The active-domain representation is stored separately by
/// the owning element type.
/// @tparam TIntegrationPoint Volume integration-point type.
/// @tparam TBoundaryIntegrationPoint Boundary integration-point type.
template<typename TIntegrationPoint, typename TBoundaryIntegrationPoint>
struct ElementCore
{
    ///@name Type definitions
    ///@{
    using IntegrationPointType = TIntegrationPoint;
    using BoundaryIntegrationPointType = TBoundaryIntegrationPoint;
    using IntegrationPointVectorType = std::vector<IntegrationPointType>;

    ///@}
    ///@name Member variables
    ///@{

    IndexType id{};
    ElementBounds bounds{};
    IntegrationPointVectorType integration_points{};
    ///@}
};
///@}
}// namespace queso
