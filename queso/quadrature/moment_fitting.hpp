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
#include <optional>
#include <vector>

//// Project includes
#include "queso/embedding/octree.h"
#include "queso/quadrature/moment_fitting_assembly.hpp"
#include "queso/quadrature/moment_fitting_types.hpp"
#include "queso/solvers/nnls.h"

namespace queso::quadrature::moment_fitting {

/// @brief Parameters for moment-fitting quadrature on a trimmed element.
struct Parameters
{
    Vector3i integration_order;  ///< Directional polynomial order.
    double residual;  ///< Target relative residual.
    std::optional<double> fictitious_domain_alpha = std::nullopt;  ///< Optional fictitious-domain weight scale.
    IndexType echo_level = 0;  ///< Diagnostic verbosity.
};

namespace detail {

    /// @brief Mutable state controlling the point-elimination loop.
    struct PointEliminationState
    {
        double residual = MIND;
        double previous_residual = 0.0;
        SizeType iteration = 0UL;
        bool point_was_eliminated = false;
    };

    /// @brief Decides whether the point-elimination loop should continue.
    /// @param rState Current point-elimination state.
    /// @param TargetResidual Requested relative residual.
    /// @param MaximumIteration Maximum number of point-elimination iterations.
    /// @return `true` if another moment-fitting/elimination iteration should be run.
    [[nodiscard]] inline bool
        ContinuePointElimination(const PointEliminationState& rState, double TargetResidual, SizeType MaximumIteration)
    {
        const bool residual_is_below_target = rState.residual < TargetResidual;
        const bool iteration_limit_is_not_reached = rState.iteration < MaximumIteration;
        return rState.point_was_eliminated || (residual_is_below_target && iteration_limit_is_not_reached);
    }

    /// @brief Erases integration points whose weight is below a tolerance.
    /// @details Uses `std::erase_if`, preserving the relative order of the remaining points.
    /// @param rPoints Integration points to modify.
    /// @param Tolerance Weight cutoff.
    /// @return `true` if at least one point was erased.
    template<typename TIntegrationPointVectorType>
    bool EraseSmallWeightPoints(TIntegrationPointVectorType& rPoints, double Tolerance = ZEROTOL)
    {
        const auto old_size = rPoints.size();
        std::erase_if(rPoints, [Tolerance](const auto& rPoint) { return rPoint.Weight() < Tolerance; });
        return rPoints.size() != old_size;
    }

    /// @brief Keeps only the largest-weight integration points.
    /// @details Uses `std::nth_element` and does not fully sort the remaining point set.
    /// @param rPoints Integration points to reduce.
    /// @param NumberOfPoints Number of largest-weight points to keep.
    /// @return `true` if points were erased.
    template<typename TIntegrationPointVectorType>
    bool KeepLargestWeightPoints(TIntegrationPointVectorType& rPoints, IndexType NumberOfPoints)
    {
        if (rPoints.size() <= NumberOfPoints) { return false; }

        const auto keep_end = rPoints.begin() + static_cast<std::ptrdiff_t>(NumberOfPoints);
        std::nth_element(rPoints.begin(), keep_end - 1, rPoints.end(), [](const auto& rPointA, const auto& rPointB) {
            return rPointA.Weight() > rPointB.Weight();
        });
        rPoints.erase(keep_end, rPoints.end());
        return true;
    }

    /// @brief Applies the first point-elimination step.
    /// @details The initial step removes near-zero weights and caps the point set to the number of moment functions.
    /// @param rPoints Integration points to reduce.
    /// @param rOrderInfo Moment-fitting order metadata.
    /// @param rWeightStats Weight statistics from the current moment-fitting solve.
    /// @return `true` if points were erased.
    template<typename TIntegrationPointVectorType>
    bool EliminateInitialPoints(
        TIntegrationPointVectorType& rPoints,
        const IntegrationOrderInfo& rOrderInfo,
        const WeightStats& rWeightStats
    )
    {
        bool points_were_eliminated = false;
        if (rWeightStats.number_of_small_weights > 0) { points_were_eliminated = EraseSmallWeightPoints(rPoints); }

        const bool largest_points_were_kept = KeepLargestWeightPoints(rPoints, rOrderInfo.number_of_functions);
        return points_were_eliminated || largest_points_were_kept;
    }

    /// @brief Eliminates low-weight points after the initial reduction step.
    /// @details Removes all points below a relative cutoff if such points exist; otherwise removes the single
    /// minimum-weight point.
    /// @param rPoints Integration points to reduce.
    /// @param rWeightStats Weight statistics from the current moment-fitting solve.
    /// @return `true` if points were erased.
    template<typename TIntegrationPointVectorType>
    bool EliminateLowWeightPoints(TIntegrationPointVectorType& rPoints, const WeightStats& rWeightStats)
    {
        if (rPoints.empty()) { return false; }

        constexpr double cut_off_ratio = 1e-8;
        const double tolerance = cut_off_ratio * rWeightStats.max_weight;
        if (rWeightStats.min_weight < tolerance) { return EraseSmallWeightPoints(rPoints, tolerance); }

        rPoints.erase(rPoints.begin() + static_cast<std::ptrdiff_t>(rWeightStats.min_weight_index));
        return true;
    }

    /// @brief Applies the point-elimination strategy for the current iteration.
    /// @param rPoints Integration points to reduce.
    /// @param rOrderInfo Moment-fitting order metadata.
    /// @param rState Current point-elimination state.
    /// @param rWeightStats Weight statistics from the current moment-fitting solve.
    /// @return `true` if points were erased.
    template<typename TIntegrationPointVectorType>
    bool EliminatePoints(
        TIntegrationPointVectorType& rPoints,
        const IntegrationOrderInfo& rOrderInfo,
        const PointEliminationState& rState,
        const WeightStats& rWeightStats
    )
    {
        if (rPoints.empty()) { return false; }

        if (rState.iteration == 0UL) { return EliminateInitialPoints(rPoints, rOrderInfo, rWeightStats); }

        return EliminateLowWeightPoints(rPoints, rWeightStats);
    }

    /// @brief Appends distributed points within trimmed domain using an octree. In each leaf node, Gauss points
    /// according to rIntegrationOrder are generated.
    ///        Only points inside the trimmed domain are considered.
    ///        Every time this function is called the octree is refined and more points are distributed.
    /// @param[out] rIntegrationPoint Existing points are preserved; new distributed points are appended.
    /// @param rOctree
    /// @param MinNumPoints Minimum number of new points to append.
    /// @param rIntegrationOrder Order of Gauss quadrature.
    template<typename TElementType, typename TOctreeOperator>
    void DistributeIntegrationPoints(
        typename TElementType::IntegrationPointVectorType& rIntegrationPoint,
        Octree<TOctreeOperator>& rOctree,
        SizeType MinNumPoints,
        const Vector3i& rIntegrationOrder
    )
    {
        constexpr IndexType max_iteration = 5;
        const SizeType target_number_of_points = rIntegrationPoint.size() + MinNumPoints;
        IndexType iteration = 0;
        while (rIntegrationPoint.size() < target_number_of_points && iteration < max_iteration) {
            const IndexType new_max_refinement = rOctree.MaxRefinementLevel() + 1;
            rOctree.Refine(std::min(new_max_refinement, IndexType{ 3 }), new_max_refinement);
            rOctree.template AddIntegrationPoints<TElementType>(rIntegrationPoint, rIntegrationOrder);
            iteration++;
        }
    }

    /// @brief Initializes immutable problem data for moment fitting.
    /// @details Retrieves active-domain boundary integration points, computes the constant terms and caches their norm.
    /// @param rElement Element whose active boundary defines the moment-fitting RHS.
    /// @param rOrderInfo Moment-fitting order metadata.
    /// @return Initialized moment-fitting problem.
    template<typename TElementType>
    MomentFittingProblem
        InitializeMomentFittingProblem(const TElementType& rElement, const IntegrationOrderInfo& rOrderInfo)
    {
        const auto boundary_ips = rElement.template GetActiveDomainBoundaryIps<
            typename TElementType::BoundaryIntegrationPointType,
            CoordinateSpace::global>();
        const auto bounds_xyz = rElement.template GetCellBounds<CoordinateSpace::global>();

        return MomentFittingProblem(ComputeConstantTerms(boundary_ips, bounds_xyz, rOrderInfo));
    }

    /// @brief Solves the moment-fitting equation and assigns non-negative weights.
    /// @details Assembles the fitting matrix, solves the NNLS problem, writes the solved weights to the integration
    /// points and collects weight statistics for point elimination.
    /// @param rIntegrationPoint Candidate integration points whose weights are overwritten.
    /// @param rGeometry Parametric geometry and Jacobian determinant for weight scaling.
    /// @param rOrderInfo Moment-fitting order metadata.
    /// @param rProblem Constant terms and cached RHS norm.
    /// @param rScratch Reusable buffers for matrix assembly and NNLS solve.
    /// @return Relative residual and solved-weight statistics.
    template<typename TElementType>
    MomentFittingResult MomentFitting(
        typename TElementType::IntegrationPointVectorType& rIntegrationPoint,
        const IntegrationGeometry& rGeometry,
        const IntegrationOrderInfo& rOrderInfo,
        const MomentFittingProblem& rProblem,
        MomentFittingScratch& rScratch
    )
    {
        const IndexType number_reduced_points = rIntegrationPoint.size();
        QuESo_ERROR_IF(number_reduced_points == 0) << "Moment fitting requires at least one integration point.\n";
        QuESo_ERROR_IF(rProblem.constant_terms_l2_norm <= ZEROTOL)
            << "Moment-fitting constant terms must have a non-zero L2 norm.\n";
        QuESo_ERROR_IF(std::abs(rGeometry.det_j) <= ZEROTOL)
            << "Moment fitting requires a non-zero Jacobian determinant.\n";

        // Assemble moment fitting matrix.
        AssembleFittingMatrix(rScratch.fitting_matrix, rIntegrationPoint, rGeometry.bounds_param, rOrderInfo);

        // Solve non-negative Least-Square-Error problem.
        rScratch.weights.resize(number_reduced_points);
        rScratch.rhs = rProblem.constant_terms;  // NNLS::solve does modify input. Therefore, copy is required.
        const double rel_residual =
            NNLS::solve(rScratch.fitting_matrix, rScratch.rhs, rScratch.weights) / rProblem.constant_terms_l2_norm;

        // Write computed weights onto integration points
        const double det_j_inv = 1.0 / rGeometry.det_j;
        WeightStats weight_stats{};
        auto weights_it = rScratch.weights.begin();
        for (IndexType point_index = 0; auto& r_point : rIntegrationPoint) {
            const double new_weight = *weights_it++ * det_j_inv;
            r_point.SetWeight(new_weight);
            weight_stats.Update(new_weight, point_index++);
        }

        return MomentFittingResult{ .residual = rel_residual, .weight_stats = weight_stats };
    }

    /// @brief Runs iterative point elimination on a candidate quadrature rule.
    /// @details Repeatedly solves moment fitting, removes low-weight points and restores the previous accepted point
    /// set if the latest reduction overshoots the target residual.
    /// @param rIntegrationPoints Candidate integration points mutated in place.
    /// @param rGeometry Parametric geometry and Jacobian determinant for weight scaling.
    /// @param rOrderInfo Moment-fitting order metadata.
    /// @param rProblem Constant terms and cached RHS norm.
    /// @param TargetResidual Target relative residual.
    /// @return Achieved relative residual.
    template<typename TElementType>
    double PointElimination(
        typename TElementType::IntegrationPointVectorType& rIntegrationPoints,
        const IntegrationGeometry& rGeometry,
        const IntegrationOrderInfo& rOrderInfo,
        const MomentFittingProblem& rProblem,
        double TargetResidual
    )
    {
        constexpr SizeType max_iteration = 10;
        PointEliminationState state;
        MomentFittingScratch scratch{};
        typename TElementType::IntegrationPointVectorType previous_integration_points{};
        previous_integration_points.reserve(rIntegrationPoints.size());

        // If any point is eliminated, we must run another moment fitting loop, to guarantee that the weights are
        // correct. Also keep iterating, until TargetResidual is stepped over.
        while (ContinuePointElimination(state, TargetResidual, max_iteration)) {
            state.point_was_eliminated = false;
            const auto fitting_result =
                MomentFitting<TElementType>(rIntegrationPoints, rGeometry, rOrderInfo, rProblem, scratch);
            state.residual = fitting_result.residual;

            if (state.iteration == 0 || state.residual < TargetResidual) {
                // Store the last accepted reduced rule before attempting another elimination.
                if (state.iteration > 0) {
                    previous_integration_points = rIntegrationPoints;
                    state.previous_residual = state.residual;
                }
                state.point_was_eliminated =
                    EliminatePoints(rIntegrationPoints, rOrderInfo, state, fitting_result.weight_stats);
            }
            state.iteration++;
        }

        if (state.point_was_eliminated && !rIntegrationPoints.empty()) {
            const auto fitting_result =
                MomentFitting<TElementType>(rIntegrationPoints, rGeometry, rOrderInfo, rProblem, scratch);
            state.residual = fitting_result.residual;
        }

        if (rIntegrationPoints.empty() && !previous_integration_points.empty()) {
            rIntegrationPoints = std::move(previous_integration_points);
            return state.previous_residual;
        }

        if (state.residual >= TargetResidual && !previous_integration_points.empty()) {
            rIntegrationPoints = std::move(previous_integration_points);
            return state.previous_residual;
        }

        return state.residual;
    }

    /// @brief Adds tensor-product quadrature points in the fictitious domain.
    /// @details Points are appended only where the parametric query point lies outside the active domain. Their weights
    /// are scaled by `Alpha`.
    /// @param rElement Element whose integration points are extended.
    /// @param rGeometry Integration geometry.
    /// @param rIntegrationOrder Order of quadrature rule.
    /// @param Alpha value that is applied to the integration weights.
    template<typename TElementType>
    void AddFictitiousIPs(
        TElementType& rElement,
        const IntegrationGeometry& rGeometry,
        const Vector3i& rIntegrationOrder,
        double Alpha
    )
    {
        auto& r_integration_points = rElement.GetIntegrationPoints();

        const auto& r_ip_list_u = IntegrationPointFactory1D::GetGauss(rIntegrationOrder[0], IntegrationMethod::gauss);
        const auto& r_ip_list_v = IntegrationPointFactory1D::GetGauss(rIntegrationOrder[1], IntegrationMethod::gauss);
        const auto& r_ip_list_w = IntegrationPointFactory1D::GetGauss(rIntegrationOrder[2], IntegrationMethod::gauss);

        const SizeType n_point = r_ip_list_u.size() * r_ip_list_v.size() * r_ip_list_w.size();
        r_integration_points.reserve(r_integration_points.size() + n_point);

        const auto& lower = rGeometry.bounds_param.lower;
        const auto& upper = rGeometry.bounds_param.upper;

        const double lower_u = lower[0];
        const double lower_v = lower[1];
        const double lower_w = lower[2];

        const double length_u = std::abs(upper[0] - lower_u);
        const double length_v = std::abs(upper[1] - lower_v);
        const double length_w = std::abs(upper[2] - lower_w);

        for (const auto& r_ip_u : r_ip_list_u) {
            const double u = lower_u + length_u * r_ip_u[0];
            const double weight_u = r_ip_u[1] * length_u;

            for (const auto& r_ip_v : r_ip_list_v) {
                const double v = lower_v + length_v * r_ip_v[0];
                const double weight_uv = weight_u * r_ip_v[1] * length_v;

                for (const auto& r_ip_w : r_ip_list_w) {
                    const PointType query_point{ u, v, lower_w + length_w * r_ip_w[0] };
                    if (!rElement.template IsInsideActiveDomain<CoordinateSpace::parametric>(query_point)) {
                        r_integration_points.emplace_back(query_point, weight_uv * r_ip_w[1] * length_w * Alpha);
                    }
                }
            }
        }
    }

}  // namespace detail

/// @brief Creates moment-fitting quadrature points for a trimmed element.
/// @details 1. Distributes initial integration points uniformly in the trimmed domain.
///          2. Computes constant terms of the moment-fitting equation.
///          3. Solves the moment-fitting equation in an iterative point-elimination algorithm.
/// See: M. Meßmer et. al: Efficient CAD-integrated isogeometric analysis of trimmed solids,
///      Comput. Methods Appl. Mech. Engrg. 400 (2022) 115584, https://doi.org/10.1016/j.cma.2022.115584.
/// @param rElement Trimmed element whose integration points are replaced.
/// @param rParameters Moment-fitting parameters.
/// @return Achieved relative moment-fitting residual, or `std::nullopt` if no valid integration points remain.
template<typename TElementType>
std::optional<double> Compute(TElementType& rElement, const Parameters& rParameters)
{
    QuESo_ERROR_IF(
        rParameters.fictitious_domain_alpha
        && (*rParameters.fictitious_domain_alpha <= 0.0 || *rParameters.fictitious_domain_alpha > 1.0)
    ) << "Fictitious-domain alpha must satisfy 0 < alpha <= 1.\n";

    const auto order_info = detail::MakeIntegrationOrderInfo(rParameters.integration_order);
    const auto problem = detail::InitializeMomentFittingProblem<TElementType>(rElement, order_info);
    const auto bounds_xyz = rElement.template GetActiveDomainBounds<CoordinateSpace::global>();
    const auto bounds_uvw = rElement.template GetActiveDomainBounds<CoordinateSpace::parametric>();
    Octree<TElementType> octree(&rElement, bounds_xyz, bounds_uvw);
    const detail::IntegrationGeometry geometry{ rElement.template GetCellBounds<CoordinateSpace::parametric>(),
                                                rElement.DetJ() };

    const SizeType max_iteration = (Math::Max(rParameters.integration_order) == 2) ? 3 : 2;
    constexpr double min_residual = 1e-2;
    double residual = MAXD;
    auto& el_integration_points = rElement.GetIntegrationPoints();
    el_integration_points.clear();
    // If residual can not be satisfied, try with more points in initial set.
    for (SizeType iteration = 0; residual > rParameters.residual && iteration < max_iteration; ++iteration) {
        // Distribute initial points via an octree.
        const SizeType min_num_points = order_info.number_of_functions * (iteration + 1);
        detail::DistributeIntegrationPoints<TElementType>(
            el_integration_points, octree, min_num_points, rParameters.integration_order
        );

        // If no point is contained in integration_points, continue.
        if (el_integration_points.empty()) { continue; }

        // Run point elimination.
        residual = detail::PointElimination<TElementType>(
            el_integration_points, geometry, order_info, problem, rParameters.residual
        );

        // If residual is very high, remove all points. Note, elements without points will be neglected.
        if (residual > min_residual) { el_integration_points.clear(); }
    }

    if (el_integration_points.empty()) {
        if (rParameters.echo_level > 2) {
            QuESo_INFO << "Warning :: Moment Fitting :: Element id: " << rElement.GetId()
                       << " neglected because no integration points remain. Residual: " << residual << ".\n";
        }
        return std::nullopt;
    }

    if (residual > rParameters.residual && rParameters.echo_level > 2) {
        QuESo_INFO << "Warning :: Moment Fitting :: Targeted residual (" << rParameters.residual
                   << ") is not achieved for element id: " << rElement.GetId() << ". Residual: " << residual << ".\n";
    }
    if (rParameters.fictitious_domain_alpha) {
        detail::AddFictitiousIPs<TElementType>(
            rElement, geometry, rParameters.integration_order, *rParameters.fictitious_domain_alpha
        );
    }
    return residual;
}

}  // namespace queso::quadrature::moment_fitting
