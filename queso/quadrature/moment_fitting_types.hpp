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
#include <functional>
#include <numeric>
#include <utility>
#include <vector>

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/solvers/nnls.h"

namespace queso::quadrature::moment_fitting::detail {

/// @brief Derived polynomial-order data used by moment-fitting assembly and point elimination.
/// @details The number of moment functions is `(order_u + 1) * (order_v + 1) * (order_w + 1)`.
struct IntegrationOrderInfo
{
    IndexType order_u;
    IndexType order_v;
    IndexType order_w;
    IndexType number_of_functions;
};

/// @brief Builds moment-fitting order metadata from the requested directional integration order.
/// @param rIntegrationOrder Polynomial/integration order in parametric u, v and w directions.
/// @return Derived integration-order metadata.
[[nodiscard]] inline IntegrationOrderInfo MakeIntegrationOrderInfo(const Vector3i& rIntegrationOrder)
{
    QuESo_ERROR_IF(rIntegrationOrder[0] < 1 || rIntegrationOrder[1] < 1 || rIntegrationOrder[2] < 1)
        << "Moment-fitting quadrature requires polynomial orders p > 0.\n";
    QuESo_ERROR_IF(rIntegrationOrder[0] > 8 || rIntegrationOrder[1] > 8 || rIntegrationOrder[2] > 8)
        << "Moment-fitting quadrature supports polynomial orders up to p=8.\n";

    return IntegrationOrderInfo{ rIntegrationOrder[0],
                                 rIntegrationOrder[1],
                                 rIntegrationOrder[2],
                                 (rIntegrationOrder[0] + 1) * (rIntegrationOrder[1] + 1) * (rIntegrationOrder[2] + 1) };
}

/// @brief Geometry data required by moment fitting in parametric space.
struct IntegrationGeometry
{
    /// Parametric bounds of the element cell.
    BoundingBoxType bounds_param;
    /// Determinant of the mapping Jacobian used to scale solved weights back to solver convention.
    double det_j;
};

/// @brief Computes the Euclidean norm of a vector.
/// @param rValues Values whose L2 norm is computed.
/// @return `sqrt(sum(value_i^2))`.
[[nodiscard]] inline double L2Norm(const std::vector<double>& rValues)
{
    return std::sqrt(std::transform_reduce(rValues.begin(), rValues.end(), 0.0, std::plus<>{}, [](double Value) {
        return Value * Value;
    }));
}

/// @brief Immutable right-hand-side data for one moment-fitting problem.
/// @details The constructor stores the constant terms and caches their L2 norm because the same RHS is used across
/// repeated NNLS solves during point elimination.
struct MomentFittingProblem
{
    /// @brief Constructs the problem from already assembled constant terms.
    /// @param ConstantTerms Moment-fitting right-hand side.
    explicit MomentFittingProblem(std::vector<double> ConstantTerms)
        : constant_terms(std::move(ConstantTerms)), constant_terms_l2_norm(L2Norm(constant_terms))
    {}

    std::vector<double> constant_terms;
    double constant_terms_l2_norm{};
};

/// @brief Reusable scratch buffers for repeated moment-fitting solves.
/// @details Reusing these buffers avoids reallocating the fitting matrix, NNLS weights and mutable RHS on every
/// point-elimination iteration.
struct MomentFittingScratch
{
    std::vector<double> weights;
    NNLS::MatrixType fitting_matrix;
    std::vector<double> rhs;
};

/// @brief Statistics of solved integration weights from one moment-fitting solve.
/// @details These values are collected while assigning weights and reused by point elimination to avoid an extra
/// min/max scan over the integration points.
struct WeightStats
{
    void Update(double Weight, IndexType Index)
    {
        max_weight = std::max(max_weight, Weight);

        if (Weight < min_weight) {
            min_weight = Weight;
            min_weight_index = Index;
        }

        if (Weight < ZEROTOL) { ++number_of_small_weights; }
    }

    double max_weight = MIND;
    double min_weight = MAXD;
    IndexType min_weight_index = 0;
    SizeType number_of_small_weights = 0;
};

/// @brief Result of one moment-fitting solve.
struct MomentFittingResult
{
    /// Relative residual `||Ax-b||_2 / ||b||_2` returned by the NNLS solve.
    double residual = 0.0;
    /// Weight statistics for the solved integration points.
    WeightStats weight_stats{};
};

}  // namespace queso::quadrature::moment_fitting::detail
