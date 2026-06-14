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
#include "queso/includes/define.hpp"
#include "queso/containers/integration_point_concepts.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  IntegrationPoint
 * @author Manuel Messmer
 * @brief Simple container that stores a 3D-point (`PointType`) and a integration weight (`double`).
 **/
class IntegrationPoint
{
public:
    ///@name Life Cycle
    ///@{

    /// 2D Constructor
    constexpr IntegrationPoint(double X, double Y, double Weight) noexcept : mPoint{ X, Y, 0.0 }, mWeight(Weight)
    {}

    /// 3D Constructor
    constexpr IntegrationPoint(double X, double Y, double Z, double Weight) noexcept
        : mPoint{ X, Y, Z }, mWeight(Weight)
    {}

    /// 3D Constructor from PointType
    constexpr IntegrationPoint(const PointType& rPoint, double Weight) noexcept : mPoint(rPoint), mWeight(Weight)
    {}

    ///@}
    ///@name Operations
    ///@{

    /// Returns underlying point coordinates.
    [[nodiscard]] constexpr const PointType& Point() const noexcept
    { return mPoint; }

    /// Access elements by index (const)
    [[nodiscard]] constexpr double operator[](IndexType i) const noexcept
    { return mPoint[i]; }

    /// Get integration weight
    [[nodiscard]] constexpr double Weight() const noexcept
    { return mWeight; }

    /// Set integration weight
    constexpr void SetWeight(double Weight) noexcept
    { mWeight = Weight; }

    ///@}
private:
    ///@name Private member variables
    ///@{
    PointType mPoint{};
    double mWeight{};
    ///@}
};// End class IntegrationPoint
///@} End QuESo classes

static_assert(concepts::IntegrationPoint<IntegrationPoint>, "IntegrationPoint must satisfy concepts::IntegrationPoint.");

}// End namespace queso
