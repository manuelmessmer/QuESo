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

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  BoundaryIntegrationPoint
 * @author Manuel Messmer
 * @brief Simple container that stores a 3D-point (`PointType`), a integration weight (`double`), and a normal vector
 * (`Vector3d`).
 **/
class BoundaryIntegrationPoint
{
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    constexpr BoundaryIntegrationPoint(double X, double Y, double Z, double Weight, const Vector3d& rNormal) noexcept
        : mPoint{ X, Y, Z }, mWeight(Weight), mNormal(rNormal)
    {}

    /// Constructor from PointType
    constexpr BoundaryIntegrationPoint(const PointType& rPoint, double Weight, const Vector3d& rNormal) noexcept
        : mPoint(rPoint), mWeight(Weight), mNormal(rNormal)
    {}

    ///@}
    ///@name Operations
    ///@{

    /// Returns underlying point coordinates.
    [[nodiscard]] constexpr const PointType& Point() const noexcept
    { return mPoint; }

    /// Access elements by index (const).
    [[nodiscard]] constexpr double operator[](IndexType i) const noexcept
    { return mPoint[i]; }

    /// Get integration weight.
    [[nodiscard]] constexpr double Weight() const noexcept
    { return mWeight; }

    /// Set integration weight.
    constexpr void SetWeight(double Weight) noexcept
    { mWeight = Weight; }

    /// @brief Get Normal vector
    /// @return
    [[nodiscard]] constexpr const Vector3d& Normal() const noexcept
    { return mNormal; }

    ///@}

private:
    ///@name Private Member Variables
    ///@{
    PointType mPoint{};
    double mWeight{};
    Vector3d mNormal{};
    ///@}
};// End class BoundaryIntegrationPoint
///@} // End QuESo classes

static_assert(std::is_trivially_copyable_v<BoundaryIntegrationPoint>,
    "BoundaryIntegrationPoint must remain trivially copyable.");

}// End namespace queso
