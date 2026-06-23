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
#include <numeric>
#include <vector>

//// Project includes
#include "queso/embedding/octree.h"
#include "queso/solvers/nnls.h"
#include "queso/utilities/polynomial_utilities.hpp"

namespace queso {

///@name QuESo Classes
///@{

/// @class  QuadratureTrimmedElement.
/// @author Manuel Messmer
/// @brief  Provides functions to create integration rules for trimmed elements.
/// @tparam TElementType
template<typename TElementType>
class QuadratureTrimmedElement
{
public:
    ///@name Type Definition
    ///@{
    using ElementType = TElementType;

    using IntegrationPointType = typename ElementType::IntegrationPointType;
    using IntegrationPointVectorType = typename ElementType::IntegrationPointVectorType;
    using IntegrationPointVectorPtrType = Unique<IntegrationPointVectorType>;

    using BoundaryIntegrationPointType = typename ElementType::BoundaryIntegrationPointType;
    using BoundaryIPsVectorType = std::vector<BoundaryIntegrationPointType>;

    using VectorType = std::vector<double>;

    ///@}
    ///@name Operations
    ///@{

    ///@brief Creates integration points for trimmed domain.
    ///@details 1. Distributes initial integration points uniformly in trimmed domain.
    ///         2. Computes constant terms of moment fitting equation.
    ///         3. Solves moment fitting equation in iterative point elimination algorithm.
    /// See: M. Meßmer et. al: Efficient CAD-integrated isogeometric analysis of trimmed solids,
    ///      Comput. Methods Appl. Mech. Engrg. 400 (2022) 115584, https://doi.org/10.1016/j.cma.2022.115584.
    ///@param rElement
    ///@param rIntegrationOrder
    ///@param Residual Targeted residual
	///@param maybeAlpha Optional alpha value for fictitious domain
    ///@param EchoLevel Default: 0
    static double AssembleIPs(
        ElementType& rElement,
        const Vector3i& rIntegrationOrder,
        double Residual,
		std::optional<double> maybeAlpha,
        IndexType EchoLevel = 0
    )
    {
        // Get boundary integration points.
        const auto boundary_ips = rElement.template GetActiveDomainBoundaryIps<
            typename TElementType::BoundaryIntegrationPointType,
            CoordinateSpace::global>();

        // Get constant terms.
        VectorType constant_terms{};
        ComputeConstantTerms(constant_terms, boundary_ips, rElement, rIntegrationOrder);

        // Construct octree. Octree is used to distribute inital points within trimmed domain.
        const auto bounds_xyz = rElement.template GetActiveDomainBounds<CoordinateSpace::global>();
        const auto bounds_uvw = rElement.template GetActiveDomainBounds<CoordinateSpace::parametric>();

        Octree<ElementType> octree(&rElement, bounds_xyz, bounds_uvw);

        // Start point elimination.
        double residual = MAXD;
        SizeType iteration = 0UL;
        SizeType point_distribution_factor = 1;
        IntegrationPointVectorType integration_points{};

        const IndexType max_iteration = (Math::Max(rIntegrationOrder) == 2) ? 4UL : 3UL;
        // If residual can not be statisfied, try with more points in initial set.
        while (residual > Residual && iteration < max_iteration) {

            // Distribute intial points via an octree.
            const SizeType min_num_points = (rIntegrationOrder[0] + 1) * (rIntegrationOrder[1] + 1)
                                            * (rIntegrationOrder[2] + 1) * (point_distribution_factor);
            DistributeIntegrationPoints(integration_points, octree, min_num_points, rIntegrationOrder);

            // If no point is contained in integration_points -> exit.
            if (integration_points.size() == 0) {
                rElement.GetIntegrationPoints().clear();
                return 1;
            }

            // Also add old, moment fitted points to new set. 'old_integration_points' only contains points with weights
            // > 0.0;
            auto& old_integration_points = rElement.GetIntegrationPoints();
            integration_points.insert(
                integration_points.end(), old_integration_points.begin(), old_integration_points.end()
            );
            old_integration_points.clear();

            // Run point elimination.
            residual = PointElimination(constant_terms, integration_points, rElement, rIntegrationOrder, Residual);

            // If residual is very high, remove all points. Note, elements without points will be neglected.
            if (residual > 1e-2) {
                auto& reduced_points = rElement.GetIntegrationPoints();
                reduced_points.clear();
            }

            // Update variables.
            point_distribution_factor *= 2;
            iteration++;
        }

        if (residual > Residual && EchoLevel > 2) {
            QuESo_INFO << "Warning :: Moment Fitting :: Targeted residual (" << Residual
                       << ") is not achieved for element id: " << rElement.GetId() << ". Residual: " << residual
                       << ".\n";
        }
        if (maybeAlpha.has_value()) { AssembleFictitiousIPs(rElement, rIntegrationOrder, maybeAlpha.value()); }
        return residual;
    }

    ///@}
protected:
    ///@name Protected Operations
    ///@{

    /// @brief Distributes point within trimmed domain using an octree. In each leaf node, Gauss points according to
    /// rIntegrationOrder are generated.
    ///        Only points inside the trimmed domain are considered.
    ///        Every time this function is called the otree is refined and more points are distributed.
    /// @param[out] rIntegrationPoint
    /// @param rOctree
    /// @param MinNumPoints Minimum Number of Points
    /// @param rIntegrationOrder Order of Gauss quadrature.
    template<typename TOctreeOperator>
    static void DistributeIntegrationPoints(
        IntegrationPointVectorType& rIntegrationPoint,
        Octree<TOctreeOperator>& rOctree,
        SizeType MinNumPoints,
        const Vector3i& rIntegrationOrder
    )
    {
        IndexType refinemen_level = rOctree.MaxRefinementLevel() + 1;
        const IndexType max_iteration = 5UL;
        IndexType iteration = 0UL;
        while (rIntegrationPoint.size() < MinNumPoints && iteration < max_iteration) {
            rOctree.Refine(std::min<IndexType>(refinemen_level, 4UL), refinemen_level);
            rIntegrationPoint.clear();
            rOctree.template AddIntegrationPoints<TElementType>(rIntegrationPoint, rIntegrationOrder);
            refinemen_level++;
            iteration++;
        }
    }

    /// @brief Computes constant terms of moment fitting equation via volume integration points.
    /// @param[out] rConstantTerms
    /// @param pIntegrationPoints (Unique<T>)
    /// @param rElement
    /// @param rIntegrationOrder
    static void ComputeConstantTerms(
        VectorType& rConstantTerms,
        const IntegrationPointVectorPtrType& pIntegrationPoints,
        const ElementType& rElement,
        const Vector3i& rIntegrationOrder
    )
    {
        // Initialize const variables.
        const PointType& a = rElement.template GetCellBounds<CoordinateSpace::parametric>().lower;
        const PointType& b = rElement.template GetCellBounds<CoordinateSpace::parametric>().upper;

        const IndexType ffactor = 1;
        const IndexType order_u = rIntegrationOrder[0];
        const IndexType order_v = rIntegrationOrder[1];
        const IndexType order_w = rIntegrationOrder[2];

        const IndexType number_of_functions =
            (order_u * ffactor + 1) * (order_v * ffactor + 1) * (order_w * ffactor + 1);

        // Resize constant terms.
        rConstantTerms.resize(number_of_functions, false);
        std::fill(rConstantTerms.begin(), rConstantTerms.end(), 0.0);

        // Loop over all boundary integration points.
        for (const auto& r_point : (*pIntegrationPoints)) {
            // For all functions
            IndexType row_index = 0UL;
            const double weight = r_point.Weight();
            for (IndexType i_x = 0; i_x <= order_u * ffactor; ++i_x) {
                for (IndexType i_y = 0; i_y <= order_v * ffactor; ++i_y) {
                    for (IndexType i_z = 0; i_z <= order_w * ffactor; ++i_z) {
                        // Assemble RHS
                        const double value = Polynomial::f_x(r_point[0], i_x, a[0], b[0])
                                             * Polynomial::f_x(r_point[1], i_y, a[1], b[1])
                                             * Polynomial::f_x(r_point[2], i_z, a[2], b[2]);
                        rConstantTerms[row_index] += value * weight;
                        row_index++;
                    }
                }
            }
        }
    }

    /// @brief Computes constant terms of moment fitting equation via boundary integration points. This functions uses
    /// the divergence theorem
    //         to transform the respective volume integrals to countour/surface integrals.
    /// @param[out] rConstantTerms
    /// @param rBoundaryIPs
    /// @param rElement
    /// @param rIntegrationOrder
    static void ComputeConstantTerms(
        VectorType& rConstantTerms,
        const BoundaryIPsVectorType& rBoundaryIPs,
        const ElementType& rElement,
        const Vector3i& rIntegrationOrder
    )
    {
        // Initialize const variables.
        const auto bounds_xyz = rElement.template GetCellBounds<CoordinateSpace::global>();

        // Constant terms / moments are evaluated in physical space.
        const PointType& a = bounds_xyz.lower;
        const PointType& b = bounds_xyz.upper;

        const IndexType ffactor = 1;
        const IndexType order_u = rIntegrationOrder[0];
        const IndexType order_v = rIntegrationOrder[1];
        const IndexType order_w = rIntegrationOrder[2];

        const IndexType number_of_functions =
            (order_u * ffactor + 1) * (order_v * ffactor + 1) * (order_w * ffactor + 1);

        // Resize constant terms.
        rConstantTerms.resize(number_of_functions, false);
        std::fill(rConstantTerms.begin(), rConstantTerms.end(), 0.0);

        /// Initialize containers for f_x and f_x_int
        // X-direction
        std::vector<double> f_x_x(order_u * ffactor + 1);
        std::vector<double> f_x_int_x(order_u * ffactor + 1);
        // Y-direction
        std::vector<double> f_x_y(order_v * ffactor + 1);
        std::vector<double> f_x_int_y(order_v * ffactor + 1);
        // Z-direction
        std::vector<double> f_x_z(order_w * ffactor + 1);
        std::vector<double> f_x_int_z(order_w * ffactor + 1);

        // Loop over all boundary integration points.
        IndexType row_index = 0;
        for (const auto& r_point : rBoundaryIPs) {
            // Note: The evaluation of polynomials is expensive. Therefore, we precompute and store values
            // for f_x_x and f_x_int at each point.
            const auto& normal = r_point.Normal();

            // X-Direction
            for (IndexType i_x = 0; i_x <= order_u * ffactor; ++i_x) {
                f_x_x[i_x] = Polynomial::f_x(r_point[0], i_x, a[0], b[0]);
                f_x_int_x[i_x] = Polynomial::f_x_int(r_point[0], i_x, a[0], b[0]);
            }
            // Y-Direction
            for (IndexType i_y = 0; i_y <= order_v * ffactor; ++i_y) {
                f_x_y[i_y] = Polynomial::f_x(r_point[1], i_y, a[1], b[1]);
                f_x_int_y[i_y] = Polynomial::f_x_int(r_point[1], i_y, a[1], b[1]);
            }
            // Z-Direction
            for (IndexType i_z = 0; i_z <= order_w * ffactor; ++i_z) {
                f_x_z[i_z] = Polynomial::f_x(r_point[2], i_z, a[2], b[2]);
                f_x_int_z[i_z] = Polynomial::f_x_int(r_point[2], i_z, a[2], b[2]);
            }

            // Assembly RHS
            row_index = 0;
            const double weight = 1.0 / 3.0 * r_point.Weight();
            for (IndexType i_x = 0; i_x <= order_u * ffactor; ++i_x) {
                for (IndexType i_y = 0; i_y <= order_v * ffactor; ++i_y) {
                    for (IndexType i_z = 0; i_z <= order_w * ffactor; ++i_z) {
                        // Compute normal for each face/triangle.
                        PointType value;
                        value[0] = f_x_int_x[i_x] * f_x_y[i_y] * f_x_z[i_z];
                        value[1] = f_x_x[i_x] * f_x_int_y[i_y] * f_x_z[i_z];
                        value[2] = f_x_x[i_x] * f_x_y[i_y] * f_x_int_z[i_z];

                        double integrand = normal[0] * value[0] + normal[1] * value[1] + normal[2] * value[2];
                        rConstantTerms[row_index] += integrand * weight;
                        row_index++;
                    }
                }
            }
        }
    }

    /// @brief Set-Up and solve moment fitting equation. Solve the moment fitting equation for the weights of the
    /// integration points.
    ///        Computed weights are directly assigned to rIntegrationPoint.
    /// @param rConstantTerms
    /// @param[out] rIntegrationPoint
    /// @param rElement
    /// @param rIntegrationOrder
    /// @return double Relative residual ||ax -b||_L2 / ||b||_L2
    static double MomentFitting(
        const VectorType& rConstantTerms,
        IntegrationPointVectorType& rIntegrationPoint,
        const ElementType& rElement,
        const Vector3i& rIntegrationOrder
    )
    {

        PointType a = rElement.template GetCellBounds<CoordinateSpace::parametric>().lower;
        PointType b = rElement.template GetCellBounds<CoordinateSpace::parametric>().upper;

        double jacobian = rElement.DetJ();

        const IndexType ffactor = 1;
        const IndexType order_u = rIntegrationOrder[0];
        const IndexType order_v = rIntegrationOrder[1];
        const IndexType order_w = rIntegrationOrder[2];

        const IndexType number_of_functions =
            (order_u * ffactor + 1) * (order_v * ffactor + 1) * (order_w * ffactor + 1);
        const IndexType number_reduced_points = rIntegrationPoint.size();

        const double l2_norm_ct =
            std::sqrt(std::inner_product(rConstantTerms.begin(), rConstantTerms.end(), rConstantTerms.begin(), 0.0));

        /// Assemble moment fitting matrix.
        NNLS::MatrixType fitting_matrix(number_of_functions * number_reduced_points, 0.0);
        const auto points_it_begin = rIntegrationPoint.begin();
        for (IndexType column_index = 0; column_index < number_reduced_points; ++column_index) {
            auto point_it = points_it_begin + static_cast<std::ptrdiff_t>(column_index);
            IndexType row_index = 0;
            for (IndexType i_x = 0; i_x <= order_u * ffactor; ++i_x) {
                for (IndexType i_y = 0; i_y <= order_v * ffactor; ++i_y) {
                    for (IndexType i_z = 0; i_z <= order_w * ffactor; ++i_z) {
                        const double value = Polynomial::f_x((*point_it)[0], i_x, a[0], b[0])
                                             * Polynomial::f_x((*point_it)[1], i_y, a[1], b[1])
                                             * Polynomial::f_x((*point_it)[2], i_z, a[2], b[2]);
                        // Matrix is serialized: Column first.
                        fitting_matrix[column_index * number_of_functions + row_index] = value;
                        row_index++;
                    }
                }
            }
        }

        // Solve non-negative Least-Square-Error problem.
        VectorType weights(number_reduced_points);
        VectorType tmp_constant_terms(rConstantTerms);// NNLS::solve does modify input. Therefore, copy is required.
        const double rel_residual = NNLS::solve(fitting_matrix, tmp_constant_terms, weights) / l2_norm_ct;

        // Write computed weights onto integration points
        for (IndexType i = 0; i < number_reduced_points; ++i) {
            // Divide by det_jacobian to account for the corresponding multiplication during the element integration
            // within the used external solver.
            double new_weight = weights[i] / (jacobian);
            rIntegrationPoint[i].SetWeight(new_weight);
        }

        return rel_residual;
    }

    /// @brief Start point elimination algorihtm. Final quadrature rule is stored in rElement.
    /// @param rConstantTerms
    /// @param rIntegrationPoint
    /// @param rElement
    /// @param rIntegrationOrder
    /// @param Residual targeted residual
    /// @return double achieved residual
    static double PointElimination(
        const VectorType& rConstantTerms,
        IntegrationPointVectorType& rIntegrationPoint,
        ElementType& rElement,
        const Vector3i& rIntegrationOrder,
        double Residual
    )
    {
        /// Initialize variables.
        const SizeType ffactor = 1;
        const SizeType order_u = rIntegrationOrder[0];
        const SizeType order_v = rIntegrationOrder[1];
        const SizeType order_w = rIntegrationOrder[2];
        const IndexType number_of_functions =
            (order_u * ffactor + 1) * (order_v * ffactor + 1) * (order_w * ffactor + 1);
        const IndexType min_number_of_points = order_u * order_v * order_w;

        const double targeted_residual = Residual;
        double global_residual = MIND;
        double prev_residual = 0.0;
        const SizeType maximum_iteration = 1000UL;
        SizeType number_iterations = 0UL;
        bool point_is_eliminated = false;
        IntegrationPointVectorType prev_solution{};

        // If any point is eliminated, we must run another moment fitting loop, to guarantee that the weights are
        // correct. Also keep iterating, until targeted_residual is stepped over.
        while (point_is_eliminated || (global_residual < targeted_residual && number_iterations < maximum_iteration)) {
            point_is_eliminated = false;
            global_residual = MomentFitting(rConstantTerms, rIntegrationPoint, rElement, rIntegrationOrder);
            if (number_iterations == 0UL) {
                /// In first iteration, revome all points but #number_of_functions
                // Sort integration points according to weight.
                std::sort(
                    rIntegrationPoint.begin(),
                    rIntegrationPoint.end(),
                    [](const IntegrationPointType& point_a, const IntegrationPointType& point_b) -> bool {
                        return point_a.Weight() > point_b.Weight();
                    }
                );
                // Only keep #number_of_functions integration points.
                rIntegrationPoint.erase(
                    rIntegrationPoint.begin() + static_cast<std::ptrdiff_t>(number_of_functions),
                    rIntegrationPoint.end()
                );

                // Additionally remove all points that are zero.
                rIntegrationPoint.erase(
                    std::remove_if(
                        rIntegrationPoint.begin(),
                        rIntegrationPoint.end(),
                        [](const IntegrationPointType& point) { return point.Weight() < ZEROTOL; }
                    ),
                    rIntegrationPoint.end()
                );

                // Stop if no points are left.
                if (rIntegrationPoint.size() == 0) break;

                point_is_eliminated = true;
            } else if (global_residual < targeted_residual && number_iterations < maximum_iteration) {
                // Store solution, in previous solution
                prev_solution.clear();
                prev_solution.insert(prev_solution.begin(), rIntegrationPoint.begin(), rIntegrationPoint.end());
                prev_residual = global_residual;

                // Find min and max weights.
                auto min_value_it = rIntegrationPoint.begin();
                double min_value = MAXD;
                double max_value = LOWESTD;
                const auto begin_it = rIntegrationPoint.begin();
                for (IndexType i = 0; i < rIntegrationPoint.size(); i++) {
                    auto it = begin_it + static_cast<std::ptrdiff_t>(i);
                    if (it->Weight() < min_value) {
                        min_value_it = it;
                        min_value = it->Weight();
                    }
                    if (it->Weight() > max_value) { max_value = it->Weight(); }
                }

                // Remove points that are very small (< EPS1*max_value)
                // However, always keep #min_number_of_points.
                const auto old_size = rIntegrationPoint.size();
                rIntegrationPoint.erase(
                    std::remove_if(
                        rIntegrationPoint.begin(),
                        rIntegrationPoint.end(),
                        [max_value](const auto& rPoint) { return rPoint.Weight() < 1e-8 * max_value; }
                    ),
                    rIntegrationPoint.end()
                );
                point_is_eliminated = (old_size != rIntegrationPoint.size());

                // If nothing was removed, remove at least one points.
                if (!point_is_eliminated && rIntegrationPoint.size() > min_number_of_points) {
                    rIntegrationPoint.erase(min_value_it);
                    point_is_eliminated = true;
                }

                // Leave loop in next iteration. Note if point_is_eliminated the moment fitting equation is solved
                // again.
                if (rIntegrationPoint.size() <= min_number_of_points) {//&& !point_is_eliminated ){
                    number_iterations = maximum_iteration + 10;
                }
            }
            number_iterations++;
        }

        auto& reduced_points = rElement.GetIntegrationPoints();
        if ((global_residual >= targeted_residual && prev_solution.size() > 0)
            || number_iterations > maximum_iteration) {
            // Return previous solution.
            reduced_points.insert(reduced_points.begin(), prev_solution.begin(), prev_solution.end());
            reduced_points.erase(
                std::remove_if(
                    reduced_points.begin(),
                    reduced_points.end(),
                    [](const IntegrationPointType& point) { return point.Weight() < ZEROTOL; }
                ),
                reduced_points.end()
            );

            return prev_residual;
        } else {
            // Return current solution.
            reduced_points.insert(reduced_points.begin(), rIntegrationPoint.begin(), rIntegrationPoint.end());
            reduced_points.erase(
                std::remove_if(
                    reduced_points.begin(),
                    reduced_points.end(),
                    [](const IntegrationPointType& point) { return point.Weight() < ZEROTOL; }
                ),
                reduced_points.end()
            );

            return global_residual;
        }
    }

    /// @brief Assembles quadrature points in the fictitious domain.
    /// @param rElement
    /// @param rIntegrationOrder Order of quadrature rule.
    /// @param Alpha value that is applied to the integration weights.
    static void AssembleFictitiousIPs(ElementType& rElement, const Vector3i& rIntegrationOrder, double Alpha)
    {
        auto& r_integration_points = rElement.GetIntegrationPoints();

        const auto& r_ip_list_u = IntegrationPointFactory1D::GetGauss(rIntegrationOrder[0], IntegrationMethod::gauss);
        const auto& r_ip_list_v = IntegrationPointFactory1D::GetGauss(rIntegrationOrder[1], IntegrationMethod::gauss);
        const auto& r_ip_list_w = IntegrationPointFactory1D::GetGauss(rIntegrationOrder[2], IntegrationMethod::gauss);

        const SizeType n_point_u = r_ip_list_u.size();
        const SizeType n_point_v = r_ip_list_v.size();
        const SizeType n_point_w = r_ip_list_w.size();

        const SizeType n_point = (n_point_u) * (n_point_v) * (n_point_w);
        r_integration_points.reserve(r_integration_points.size() + n_point);

        const PointType lower_bound_param = rElement.template GetCellBounds<CoordinateSpace::parametric>().lower;
        const PointType upper_bound_param = rElement.template GetCellBounds<CoordinateSpace::parametric>().upper;

        const double length_u = upper_bound_param[0] - lower_bound_param[0];
        const double length_v = upper_bound_param[1] - lower_bound_param[1];
        const double length_w = upper_bound_param[2] - lower_bound_param[2];
        const double scaled_volume = length_u * length_v * length_w * Alpha;

        for (const auto& ip_u : r_ip_list_u) {
            const double u_coord = lower_bound_param[0] + length_u * ip_u[0];
            for (const auto& ip_v : r_ip_list_v) {
                const double v_coord = lower_bound_param[1] + length_v * ip_v[0];
                const double uv_weight = ip_u[1] * ip_v[1];
                for (const auto& ip_w : r_ip_list_w) {
                    const PointType query_point{ u_coord, v_coord, lower_bound_param[2] + length_w * ip_w[0] };
                    if (!rElement.template IsInsideActiveDomain<CoordinateSpace::parametric>(query_point)) {
                        r_integration_points.emplace_back(query_point, uv_weight * ip_w[1] * scaled_volume);
                    }
                }
            }
        }
    }
};// End Class
}// End namespace queso
