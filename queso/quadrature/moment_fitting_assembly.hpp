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
#include <utility>
#include <vector>

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/quadrature/moment_fitting_types.hpp"
#include "queso/solvers/nnls.h"
#include "queso/utilities/polynomial_utilities.hpp"

namespace queso {
class BoundaryIntegrationPoint;
class IntegrationPoint;
}  // namespace queso

namespace queso::quadrature::moment_fitting::detail {

/// @brief Invokes a callable for each compile-time index in `[0, Count)`.
/// @param rFunction Callable receiving `std::integral_constant<IndexType, I>`.
template<IndexType Count, typename TFunction>
void StaticFor(TFunction&& rFunction)
{
    [&]<IndexType... Indices>(std::integer_sequence<IndexType, Indices...>) {
        (rFunction(std::integral_constant<IndexType, Indices>{}), ...);
    }(std::make_integer_sequence<IndexType, Count>{});
}

/// @brief Dispatches a runtime order in `[0, 8]` to a compile-time integral constant.
/// @details Callers must handle unsupported values before calling this function. The default branch is unreachable for
/// guarded calls and aborts if the precondition is violated.
/// @param Value Runtime polynomial order.
/// @param rFunction Callable receiving `std::integral_constant<IndexType, Value>`.
/// @return The callable return value.
template<typename TFunction>
decltype(auto) DispatchOrder(IndexType Value, TFunction&& rFunction)
{
    switch (Value) {
    case 0:
        return rFunction(std::integral_constant<IndexType, 0>{});
    case 1:
        return rFunction(std::integral_constant<IndexType, 1>{});
    case 2:
        return rFunction(std::integral_constant<IndexType, 2>{});
    case 3:
        return rFunction(std::integral_constant<IndexType, 3>{});
    case 4:
        return rFunction(std::integral_constant<IndexType, 4>{});
    case 5:
        return rFunction(std::integral_constant<IndexType, 5>{});
    case 6:
        return rFunction(std::integral_constant<IndexType, 6>{});
    case 7:
        return rFunction(std::integral_constant<IndexType, 7>{});
    case 8:
        return rFunction(std::integral_constant<IndexType, 8>{});
    default:
        Unreachable("Moment-fitting polynomial order dispatch received an unsupported order.");
    }
}

/// @brief Fills Legendre basis values for one coordinate direction.
/// @param Coordinate Coordinate where the basis is evaluated.
/// @param Lower Lower bound of the coordinate interval.
/// @param Upper Upper bound of the coordinate interval.
/// @param rBasis Output basis values for orders `0..Order`.
template<IndexType Order>
void FillBasisValues(double Coordinate, double Lower, double Upper, std::array<double, Order + 1>& rBasis)
{
    StaticFor<Order + 1>([&](auto I) {
        constexpr IndexType i = I;
        rBasis[i] = polynomial::f_x<i>(Coordinate, Lower, Upper);
    });
}

/// @brief Fills Legendre basis values and their integrals for one coordinate direction.
/// @param Coordinate Coordinate where the basis is evaluated.
/// @param Lower Lower bound of the coordinate interval.
/// @param Upper Upper bound of the coordinate interval.
/// @param rBasis Output basis values for orders `0..Order`.
/// @param rIntegratedBasis Output integrated basis values for orders `0..Order`.
template<IndexType Order>
void FillBasis(
    double Coordinate,
    double Lower,
    double Upper,
    std::array<double, Order + 1>& rBasis,
    std::array<double, Order + 1>& rIntegratedBasis
)
{
    StaticFor<Order + 1>([&](auto I) {
        constexpr IndexType i = I;
        rBasis[i] = polynomial::f_x<i>(Coordinate, Lower, Upper);
        rIntegratedBasis[i] = polynomial::f_x_int<i>(Coordinate, Lower, Upper);
    });
}

/// @brief Computes constant terms with compile-time polynomial orders.
/// @details Optimized path for supported orders. Uses stack-allocated basis arrays and templated polynomial evaluation.
/// @param rBoundaryIPs Boundary integration points in global coordinates.
/// @param rBoundsGlobal Global element-cell bounds used to evaluate physical-space basis functions.
/// @return Moment-fitting right-hand side.
template<typename TBoundaryIPsVectorType, IndexType OrderU, IndexType OrderV, IndexType OrderW>
std::vector<double>
    ComputeConstantTermsFixedOrder(const TBoundaryIPsVectorType& rBoundaryIPs, const BoundingBoxType& rBoundsGlobal)
{
    const PointType& a = rBoundsGlobal.lower;
    const PointType& b = rBoundsGlobal.upper;

    constexpr SizeType number_basis_x = OrderU + 1;
    constexpr SizeType number_basis_y = OrderV + 1;
    constexpr SizeType number_basis_z = OrderW + 1;
    constexpr SizeType number_of_functions = number_basis_x * number_basis_y * number_basis_z;

    std::vector<double> constant_terms(number_of_functions, 0.0);

    std::array<double, number_basis_x> basis_x{};
    std::array<double, number_basis_x> integrated_basis_x{};
    std::array<double, number_basis_y> basis_y{};
    std::array<double, number_basis_y> integrated_basis_y{};
    std::array<double, number_basis_z> basis_z{};
    std::array<double, number_basis_z> integrated_basis_z{};

    constexpr double one_third = 1.0 / 3.0;

    for (const auto& r_point : rBoundaryIPs) {
        const auto& normal = r_point.Normal();

        FillBasis<OrderU>(r_point[0], a[0], b[0], basis_x, integrated_basis_x);
        FillBasis<OrderV>(r_point[1], a[1], b[1], basis_y, integrated_basis_y);
        FillBasis<OrderW>(r_point[2], a[2], b[2], basis_z, integrated_basis_z);

        auto constant_terms_it = constant_terms.begin();
        const double weighted_area = one_third * r_point.Weight();
        const double weighted_normal_x = weighted_area * normal[0];
        const double weighted_normal_y = weighted_area * normal[1];
        const double weighted_normal_z = weighted_area * normal[2];
        for (IndexType i_x = 0; i_x < number_basis_x; ++i_x) {
            const double basis_value_x = basis_x[i_x];
            const double integrated_basis_value_x = integrated_basis_x[i_x];
            for (IndexType i_y = 0; i_y < number_basis_y; ++i_y) {
                const double basis_value_y = basis_y[i_y];
                const double integrated_basis_value_y = integrated_basis_y[i_y];
                const double z_basis_factor = weighted_normal_x * integrated_basis_value_x * basis_value_y
                                              + weighted_normal_y * basis_value_x * integrated_basis_value_y;
                const double z_integrated_basis_factor = weighted_normal_z * basis_value_x * basis_value_y;
                auto integrated_basis_z_it = integrated_basis_z.begin();
                for (const double basis_value_z : basis_z) {
                    *constant_terms_it++ +=
                        z_basis_factor * basis_value_z + z_integrated_basis_factor * *integrated_basis_z_it++;
                }
            }
        }
    }

    return constant_terms;
}

/// @brief Computes constant terms for a moment-fitting problem.
/// @details Supports orders `0..8` and dispatches to compile-time polynomial evaluation.
/// @param rBoundaryIPs Boundary integration points in global coordinates.
/// @param rBoundsGlobal Global element-cell bounds used to evaluate physical-space basis functions.
/// @param rOrderInfo Moment-fitting order metadata.
/// @return Moment-fitting right-hand side.
template<typename TBoundaryIPsVectorType>
std::vector<double> ComputeConstantTerms(
    const TBoundaryIPsVectorType& rBoundaryIPs,
    const BoundingBoxType& rBoundsGlobal,
    const IntegrationOrderInfo& rOrderInfo
)
{
    QuESo_ERROR_IF(rOrderInfo.order_u > 8 || rOrderInfo.order_v > 8 || rOrderInfo.order_w > 8)
        << "Moment-fitting quadrature supports polynomial orders up to p=8.\n";

    return DispatchOrder(rOrderInfo.order_u, [&](auto U) {
        return DispatchOrder(rOrderInfo.order_v, [&](auto V) {
            return DispatchOrder(rOrderInfo.order_w, [&](auto W) {
                constexpr IndexType order_u = decltype(U)::value;
                constexpr IndexType order_v = decltype(V)::value;
                constexpr IndexType order_w = decltype(W)::value;
                return ComputeConstantTermsFixedOrder<TBoundaryIPsVectorType, order_u, order_v, order_w>(
                    rBoundaryIPs, rBoundsGlobal
                );
            });
        });
    });
}

/// @brief Assembles the NNLS fitting matrix with compile-time polynomial orders.
/// @details Optimized path for supported orders. Uses stack-allocated basis arrays and templated polynomial evaluation.
/// The matrix is serialized column-first as expected by `NNLS::solve`.
/// @param[out] rFittingMatrix Serialized fitting matrix.
/// @param rIntegrationPoints Candidate integration points in parametric coordinates.
/// @param rBoundsParam Parametric element-cell bounds.
template<typename TIntegrationPointVectorType, IndexType OrderU, IndexType OrderV, IndexType OrderW>
void AssembleFittingMatrixFixedOrder(
    NNLS::MatrixType& rFittingMatrix,
    const TIntegrationPointVectorType& rIntegrationPoints,
    const BoundingBoxType& rBoundsParam
)
{
    const auto& a = rBoundsParam.lower;
    const auto& b = rBoundsParam.upper;

    constexpr SizeType number_basis_u = OrderU + 1;
    constexpr SizeType number_basis_v = OrderV + 1;
    constexpr SizeType number_basis_w = OrderW + 1;
    constexpr SizeType number_of_functions = number_basis_u * number_basis_v * number_basis_w;

    rFittingMatrix.resize(number_of_functions * rIntegrationPoints.size());
    std::array<double, number_basis_u> basis_u{};
    std::array<double, number_basis_v> basis_v{};
    std::array<double, number_basis_w> basis_w{};

    auto fitting_matrix_it = rFittingMatrix.begin();
    for (const auto& r_point : rIntegrationPoints) {
        FillBasisValues<OrderU>(r_point[0], a[0], b[0], basis_u);
        FillBasisValues<OrderV>(r_point[1], a[1], b[1], basis_v);
        FillBasisValues<OrderW>(r_point[2], a[2], b[2], basis_w);

        // Matrix is serialized column first.
        for (const double basis_value_u : basis_u) {
            for (const double basis_value_v : basis_v) {
                for (const double basis_value_w : basis_w) {
                    *fitting_matrix_it++ = basis_value_u * basis_value_v * basis_value_w;
                }
            }
        }
    }
}

/// @brief Assembles the NNLS fitting matrix for a moment-fitting solve.
/// @details Supports orders `0..8` and dispatches to compile-time polynomial evaluation.
/// @param[out] rFittingMatrix Serialized fitting matrix.
/// @param rIntegrationPoints Candidate integration points in parametric coordinates.
/// @param rBoundsParam Parametric element-cell bounds.
/// @param rOrderInfo Moment-fitting order metadata.
template<typename TIntegrationPointVectorType>
void AssembleFittingMatrix(
    NNLS::MatrixType& rFittingMatrix,
    const TIntegrationPointVectorType& rIntegrationPoints,
    const BoundingBoxType& rBoundsParam,
    const IntegrationOrderInfo& rOrderInfo
)
{
    QuESo_ERROR_IF(rOrderInfo.order_u > 8 || rOrderInfo.order_v > 8 || rOrderInfo.order_w > 8)
        << "Moment-fitting quadrature supports polynomial orders up to p=8.\n";

    DispatchOrder(rOrderInfo.order_u, [&](auto U) {
        DispatchOrder(rOrderInfo.order_v, [&](auto V) {
            DispatchOrder(rOrderInfo.order_w, [&](auto W) {
                constexpr IndexType order_u = decltype(U)::value;
                constexpr IndexType order_v = decltype(V)::value;
                constexpr IndexType order_w = decltype(W)::value;
                AssembleFittingMatrixFixedOrder<TIntegrationPointVectorType, order_u, order_v, order_w>(
                    rFittingMatrix, rIntegrationPoints, rBoundsParam
                );
            });
        });
    });
}

extern template std::vector<double> ComputeConstantTerms<std::vector<BoundaryIntegrationPoint>>(
    const std::vector<BoundaryIntegrationPoint>& rBoundaryIPs,
    const BoundingBoxType& rBoundsGlobal,
    const IntegrationOrderInfo& rOrderInfo
);

extern template void AssembleFittingMatrix<std::vector<IntegrationPoint>>(
    NNLS::MatrixType& rFittingMatrix,
    const std::vector<IntegrationPoint>& rIntegrationPoints,
    const BoundingBoxType& rBoundsParam,
    const IntegrationOrderInfo& rOrderInfo
);

}  // namespace queso::quadrature::moment_fitting::detail
