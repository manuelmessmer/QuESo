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

/// @brief Provides functions to evaluate Legendre polynomials and their integrals.
namespace queso::polynomial {

/// @brief Computes `x` raised to a compile-time power.
template<std::size_t Power>
[[nodiscard]] inline double power(double x) noexcept
{
    if constexpr (Power == 0) {
        return 1.0;
    } else if constexpr (Power == 1) {
        return x;
    } else {
        return x * power<Power - 1>(x);
    }
}

/// @brief Evaluates the Legendre basis function of order `Order` on interval `(a,b)`.
template<IndexType Order>
[[nodiscard]] inline double f_x(double x, double a, double b) noexcept
{
    const double tmp_x = (2.0 * x - a - b) / (b - a);
    if constexpr (Order == 0) {
        return 1.0;
    } else if constexpr (Order == 1) {
        return tmp_x;
    } else if constexpr (Order == 2) {
        return 1.0 / 2.0 * (3.0 * power<2>(tmp_x) - 1.0);
    } else if constexpr (Order == 3) {
        return 1.0 / 2.0 * (5.0 * power<3>(tmp_x) - 3.0 * tmp_x);
    } else if constexpr (Order == 4) {
        return 1.0 / 8.0 * (35.0 * power<4>(tmp_x) - 30.0 * power<2>(tmp_x) + 3.0);
    } else if constexpr (Order == 5) {
        return 1.0 / 8.0 * (63.0 * power<5>(tmp_x) - 70.0 * power<3>(tmp_x) + 15.0 * tmp_x);
    } else if constexpr (Order == 6) {
        return 1.0 / 16.0 * (231.0 * power<6>(tmp_x) - 315.0 * power<4>(tmp_x) + 105.0 * power<2>(tmp_x) - 5.0);
    } else if constexpr (Order == 7) {
        return 1.0 / 16.0
               * (429.0 * power<7>(tmp_x) - 693.0 * power<5>(tmp_x) + 315.0 * power<3>(tmp_x) - 35.0 * tmp_x);
    } else if constexpr (Order == 8) {
        return 1.0 / 128.0
               * (6435.0 * power<8>(tmp_x) - 12012.0 * power<6>(tmp_x) + 6930.0 * power<4>(tmp_x)
                  - 1260.0 * power<2>(tmp_x) + 35.0);
    } else {
        static_assert(always_false_v<std::integral_constant<IndexType, Order>>, "Supported Legendre orders are 0..8.");
    }
}

/// @brief Evaluates the antiderivative of the Legendre basis function of order `Order` on interval `(a,b)`.
template<IndexType Order>
[[nodiscard]] inline double f_x_int(double x, double a, double b) noexcept
{
    const double shifted_x = a + b - 2.0 * x;
    const double length = a - b;
    if constexpr (Order == 0) {
        return x;
    } else if constexpr (Order == 1) {
        return -power<2>(shifted_x) / (4.0 * length);
    } else if constexpr (Order == 2) {
        return -x / 2.0 - power<3>(shifted_x) / (4.0 * power<2>(length));
    } else if constexpr (Order == 3) {
        return (3.0 * power<2>(shifted_x)) / (8.0 * length) - (5.0 * power<4>(shifted_x)) / (16.0 * power<3>(length));
    } else if constexpr (Order == 4) {
        return (3.0 * x) / 8.0 + (5.0 * power<3>(shifted_x)) / (8.0 * power<2>(length))
               - (7.0 * power<5>(shifted_x)) / (16.0 * power<4>(length));
    } else if constexpr (Order == 5) {
        return (35.0 * power<4>(shifted_x)) / (32.0 * power<3>(length)) - (15.0 * power<2>(shifted_x)) / (32.0 * length)
               - (21.0 * power<6>(shifted_x)) / (32.0 * power<5>(length));
    } else if constexpr (Order == 6) {
        return (63.0 * power<5>(shifted_x)) / (32.0 * power<4>(length))
               - (35.0 * power<3>(shifted_x)) / (32.0 * power<2>(length)) - (5.0 * x) / 16.0
               - (33.0 * power<7>(shifted_x)) / (32.0 * power<6>(length));
    } else if constexpr (Order == 7) {
        return (35.0 * power<2>(shifted_x)) / (64.0 * length)
               - (315.0 * power<4>(shifted_x)) / (128.0 * power<3>(length))
               + (231.0 * power<6>(shifted_x)) / (64.0 * power<5>(length))
               - (429.0 * power<8>(shifted_x)) / (256.0 * power<7>(length));
    } else if constexpr (Order == 8) {
        return (35.0 * x) / 128.0 + (105.0 * power<3>(shifted_x)) / (64.0 * power<2>(length))
               - (693.0 * power<5>(shifted_x)) / (128.0 * power<4>(length))
               + (429.0 * power<7>(shifted_x)) / (64.0 * power<6>(length))
               - (715.0 * power<9>(shifted_x)) / (256.0 * power<8>(length));
    } else {
        static_assert(always_false_v<std::integral_constant<IndexType, Order>>, "Supported Legendre orders are 0..8.");
    }
}

}  // namespace queso::polynomial
