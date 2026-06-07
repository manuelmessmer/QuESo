//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \\'
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
#include <concepts>
#include <type_traits>

//// Project includes
#include "queso/includes/define.hpp"

namespace queso::concepts {

template<typename T>
concept IntegrationPoint =
    std::is_trivially_copyable_v<T> && requires(T point, const T const_point, IndexType i, double weight) {
        { const_point.Point() } -> std::same_as<const PointType&>;
        { const_point[i] } -> std::same_as<double>;
        { const_point.Weight() } -> std::same_as<double>;
        { point.SetWeight(weight) } -> std::same_as<void>;
    };

template<typename T>
concept BoundaryIntegrationPoint = IntegrationPoint<T> && requires(const T const_point) {
    { const_point.Normal() } -> std::same_as<const Vector3d&>;
};

}// namespace queso::concepts
