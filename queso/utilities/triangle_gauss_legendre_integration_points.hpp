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
#include <array>
#include <span>

//// Project includes
#include "queso/containers/integration_point.hpp"
#include "queso/includes/exception.hpp"

namespace queso::triangle_utilities {

/// Selects the maximum polynomial degree integrated exactly by triangle quadrature.
enum class TriangleQuadratureOrder {
    Degree1 = 1,  ///< One-point rule, exact through degree 1.
    Degree2 = 2,  ///< Three-point rule, exact through degree 2.
    Degree4 = 4,  ///< Six-point rule, exact through degree 4.
    Degree6 = 6  ///< Twelve-point rule, exact through degree 6.
};

namespace detail {

    /// Returns linear shape-function values at reference-triangle coordinates (xi, eta, 0).
    [[nodiscard]] constexpr std::array<double, 3> ShapeFunctionValues(const PointType& rPoint) noexcept
    { return { 1.0 - rPoint[0] - rPoint[1], rPoint[0], rPoint[1] }; }

    /// Stores a rule and its matching shape functions, evaluated together at compile time.
    template<SizeType TSize>
    struct ReferenceRule
    {
        std::array<IntegrationPoint, TSize> points;
        std::array<std::array<double, 3>, TSize> shape_functions;

        constexpr explicit ReferenceRule(const std::array<IntegrationPoint, TSize>& rPoints) noexcept
            : points(rPoints), shape_functions{}
        {
            for (SizeType I = 0; I < TSize; ++I) { shape_functions[I] = ShapeFunctionValues(points[I].Point()); }
        }
    };

    inline constexpr ReferenceRule Rule1{ std::array{ IntegrationPoint{ 1.0 / 3.0, 1.0 / 3.0, 0.5 } } };

    inline constexpr ReferenceRule Rule2{ std::array{ IntegrationPoint{ 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0 },
                                                      IntegrationPoint{ 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0 },
                                                      IntegrationPoint{ 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 } } };

    inline constexpr ReferenceRule Rule3 = [] {
        constexpr double Wa = 0.054975871827661;
        constexpr double Wb = 0.1116907948390055;
        constexpr double Na1 = 0.816847572980459;
        constexpr double Nb1 = 0.108103018168070;
        constexpr double Na2 = 0.091576213509771;
        constexpr double Nb2 = 0.445948490915965;
        return ReferenceRule{ std::array{ IntegrationPoint{ Na2, Na2, Wa },
                                          IntegrationPoint{ Na1, Na2, Wa },
                                          IntegrationPoint{ Na2, Na1, Wa },
                                          IntegrationPoint{ Nb2, Nb2, Wb },
                                          IntegrationPoint{ Nb1, Nb2, Wb },
                                          IntegrationPoint{ Nb2, Nb1, Wb } } };
    }();

    inline constexpr ReferenceRule Rule4 = [] {
        constexpr double Wa = 0.025422453185103408460;
        constexpr double Wb = 0.058393137863189683013;
        constexpr double Wc = 0.041425537809186787597;
        constexpr double N1 = 0.87382197101699554332;
        constexpr double N2 = 0.063089014491502228340;
        constexpr double N3 = 0.50142650965817915742;
        constexpr double N4 = 0.24928674517091042129;
        constexpr double N5 = 0.053145049844816947353;
        constexpr double N6 = 0.31035245103378440542;
        constexpr double N7 = 0.63650249912139864723;
        return ReferenceRule{ std::array{ IntegrationPoint{ N1, N2, Wa },
                                          IntegrationPoint{ N2, N1, Wa },
                                          IntegrationPoint{ N2, N2, Wa },
                                          IntegrationPoint{ N3, N4, Wb },
                                          IntegrationPoint{ N4, N3, Wb },
                                          IntegrationPoint{ N4, N4, Wb },
                                          IntegrationPoint{ N5, N6, Wc },
                                          IntegrationPoint{ N6, N5, Wc },
                                          IntegrationPoint{ N5, N7, Wc },
                                          IntegrationPoint{ N6, N7, Wc },
                                          IntegrationPoint{ N7, N5, Wc },
                                          IntegrationPoint{ N7, N6, Wc } } };
    }();

    /// Non-owning views into a rule's static point and precomputed shape-function tables.
    struct ReferenceRuleView
    {
        std::span<const IntegrationPoint> points;
        std::span<const std::array<double, 3>> shape_functions;
    };

    [[nodiscard]] constexpr ReferenceRuleView GetRule(TriangleQuadratureOrder Order)
    {
        switch (Order) {
        case TriangleQuadratureOrder::Degree1:
            return { Rule1.points, Rule1.shape_functions };
        case TriangleQuadratureOrder::Degree2:
            return { Rule2.points, Rule2.shape_functions };
        case TriangleQuadratureOrder::Degree4:
            return { Rule3.points, Rule3.shape_functions };
        case TriangleQuadratureOrder::Degree6:
            return { Rule4.points, Rule4.shape_functions };
        }
        QuESo_ERROR << "Unsupported triangle quadrature order: " << static_cast<int>(Order)
                    << ". Expected Degree1, Degree2, Degree4, or Degree6.";
    }

}  // namespace detail

/// @brief Returns quadrature points on the reference triangle (0, 0), (1, 0), (0, 1).
/// @details Coordinates are (xi, eta, 0). Weights sum to 1/2. Orders Degree1, Degree2, Degree4, and Degree6
/// use 1, 3, 6, and 12 points, respectively. Valid lookups support constant evaluation.
/// @param Order Maximum polynomial degree integrated exactly.
/// @return Read-only span into immutable static storage, valid for the program lifetime.
/// @throws Exception If Order is not a supported enumerator; invalid orders cannot be used in constant expressions.
[[nodiscard]] constexpr std::span<const IntegrationPoint> GetIntegrationPoints(TriangleQuadratureOrder Order)
{ return detail::GetRule(Order).points; }

/// @brief Returns the number of points in the selected triangle quadrature rule.
/// @param Order Maximum polynomial degree integrated exactly.
/// @return Number of integration points; supports constant evaluation.
/// @throws Exception If Order is not a supported enumerator.
[[nodiscard]] constexpr SizeType NumberOfIntegrationPoints(TriangleQuadratureOrder Order)
{ return GetIntegrationPoints(Order).size(); }

}  // namespace queso::triangle_utilities
