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
#include <algorithm>
#include <cmath>
#include <limits>

//// Project includes
#include "queso/includes/core_definitions.hpp"
#include "queso/includes/exception.hpp"

namespace queso {

/// @brief Immutable model-scale tolerances for geometric lengths, areas, and volumes.
/// @details High-level owners retain this coordinated policy. Low-level predicates receive resolved scalar values when
/// a single threshold is their complete physical contract.
class GeometryTolerance
{
public:
    /// @brief Inputs used to resolve model-scale geometric tolerances.
    struct Scale
    {
        double length_scale;  ///< Characteristic geometric length.
        double coordinate_scale;  ///< Maximum absolute coordinate magnitude.
    };

    /// @brief Creates tolerances resolved from geometric and coordinate scales.
    /// @param ScaleInput Named geometric and coordinate scales.
    /// @return Resolved geometry tolerance.
    [[nodiscard]] static GeometryTolerance FromScale(Scale ScaleInput)
    {
        QuESo_ERROR_IF(!std::isfinite(ScaleInput.length_scale) || ScaleInput.length_scale <= 0.0)
            << "Geometry-tolerance length scale must be finite and positive.\n";
        QuESo_ERROR_IF(!std::isfinite(ScaleInput.coordinate_scale) || ScaleInput.coordinate_scale <= 0.0)
            << "Geometry-tolerance coordinate scale must be finite and positive.\n";

        const double coordinate_scale = std::max(1.0, ScaleInput.coordinate_scale);
        const double coordinate_resolution =
            CoordinateResolutionMultiplier * std::numeric_limits<double>::epsilon() * coordinate_scale;
        const double snap_distance = std::max(SnapRelativeCoefficient * ScaleInput.length_scale, coordinate_resolution);
        const double zero_length = std::max(ZeroRelativeCoefficient * ScaleInput.length_scale, coordinate_resolution);
        const double zero_area = zero_length * zero_length;
        const double zero_volume = zero_area * zero_length;

        QuESo_ERROR_IF(
            !std::isfinite(snap_distance) || !std::isfinite(zero_length) || !std::isfinite(zero_area)
            || !std::isfinite(zero_volume) || snap_distance >= ScaleInput.length_scale
            || zero_length >= ScaleInput.length_scale || zero_area >= ScaleInput.length_scale * ScaleInput.length_scale
            || zero_volume >= ScaleInput.length_scale * ScaleInput.length_scale * ScaleInput.length_scale
        ) << "Geometry tolerance cannot be represented safely at the supplied scales.\n";

        return GeometryTolerance(snap_distance, zero_length, zero_area, zero_volume);
    }

    /// @brief Creates tolerances resolved from axis-aligned bounds.
    /// @param rBounds Finite bounds with strictly positive extent along every axis.
    /// @return Resolved geometry tolerance for the bounds' coordinate space.
    [[nodiscard]] static GeometryTolerance FromBounds(const BoundingBoxType& rBounds)
    {
        PointType delta{};
        for (IndexType axis = 0; axis < 3; ++axis) {
            QuESo_ERROR_IF(
                !std::isfinite(rBounds.lower[axis]) || !std::isfinite(rBounds.upper[axis])
                || rBounds.lower[axis] >= rBounds.upper[axis]
            ) << "Geometry-tolerance bounds must be finite with positive extent.\n";
            delta[axis] = rBounds.upper[axis] - rBounds.lower[axis];
        }

        return FromScale(
            { .length_scale = std::max({ delta[0], delta[1], delta[2] }),
              .coordinate_scale = std::max(
                  { 1.0,
                    std::abs(rBounds.lower[0]),
                    std::abs(rBounds.lower[1]),
                    std::abs(rBounds.lower[2]),
                    std::abs(rBounds.upper[0]),
                    std::abs(rBounds.upper[1]),
                    std::abs(rBounds.upper[2]) }
              ) }
        );
    }

    /// @brief Returns the coordinate snapping and point-welding distance.
    /// @return Absolute length in global coordinates.
    [[nodiscard]] double SnapDistance() const noexcept
    { return mSnapDistance; }

    /// @brief Returns the geometric zero-length threshold.
    /// @return Absolute length in global coordinates.
    [[nodiscard]] double ZeroLength() const noexcept
    { return mZeroLength; }

    /// @brief Returns the geometric zero-area threshold.
    /// @return Absolute area in global coordinates.
    [[nodiscard]] double ZeroArea() const noexcept
    { return mZeroArea; }

    /// @brief Returns the geometric zero-volume threshold.
    /// @return Absolute volume in global coordinates.
    [[nodiscard]] double ZeroVolume() const noexcept
    { return mZeroVolume; }

    /// @brief Tests whether two coordinates are within the snap distance.
    /// @param Left First coordinate.
    /// @param Right Second coordinate.
    /// @return True when the coordinates are the same under the snapping policy.
    [[nodiscard]] bool CoordinatesAreSame(double Left, double Right) const noexcept
    { return std::abs(Left - Right) <= mSnapDistance; }

    /// @brief Tests whether two points are within the Euclidean snap distance.
    /// @param rLeft First point.
    /// @param rRight Second point.
    /// @return True when the points are the same under the welding policy.
    [[nodiscard]] bool PointsAreSame(PointView rLeft, PointView rRight) const noexcept
    { return std::hypot(rLeft[0] - rRight[0], rLeft[1] - rRight[1], rLeft[2] - rRight[2]) <= mSnapDistance; }

    /// @brief Tests whether a signed length is geometrically zero.
    /// @param Value Signed length.
    /// @return True when the magnitude does not exceed the zero-length threshold.
    [[nodiscard]] bool IsZeroLength(double Value) const noexcept
    { return std::abs(Value) <= mZeroLength; }

    /// @brief Tests whether a signed area is geometrically zero.
    /// @param Value Signed area.
    /// @return True when the magnitude does not exceed the zero-area threshold.
    [[nodiscard]] bool IsZeroArea(double Value) const noexcept
    { return std::abs(Value) <= mZeroArea; }

    /// @brief Tests whether a signed volume is geometrically zero.
    /// @param Value Signed volume.
    /// @return True when the magnitude does not exceed the zero-volume threshold.
    [[nodiscard]] bool IsZeroVolume(double Value) const noexcept
    { return std::abs(Value) <= mZeroVolume; }

private:
    inline static constexpr double SnapRelativeCoefficient = 1e-12;
    inline static constexpr double ZeroRelativeCoefficient = 1e-14;
    inline static constexpr double CoordinateResolutionMultiplier = 32.0;

    constexpr GeometryTolerance(double SnapDistance, double ZeroLength, double ZeroArea, double ZeroVolume) noexcept
        : mSnapDistance(SnapDistance), mZeroLength(ZeroLength), mZeroArea(ZeroArea), mZeroVolume(ZeroVolume)
    {}

    double mSnapDistance;
    double mZeroLength;
    double mZeroArea;
    double mZeroVolume;
};

}  // namespace queso
