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

//// External includes
#include <doctest/doctest.h>

//// STL includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/quadrature/moment_fitting.hpp"
#include "queso/tests/cpp/helpers/trimmed_element_test_helpers.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/mesh_utilities.h"

// This suite tests public moment-fitting Compute behavior on trimmed elements. Related coverage:
// test_element_trimmed.cpp checks the TrimmedElement container API, test_moment_fitting.cpp checks NNLS solve behavior,
// and test_point_elimination.cpp checks reduced rules on larger trimmed geometries.

namespace queso::testing {

namespace {

    using namespace TrimmedElementTestHelpers;
    namespace mf = quadrature::moment_fitting;

    SizeType CountOutsideActiveDomain(const TrimmedElementType& rElement)
    {
        SizeType result = 0;
        for (const auto& rPoint : rElement.GetIntegrationPoints()) {
            if (!rElement.template IsInsideActiveDomain<CoordinateSpace::parametric>(rPoint.Point())) { ++result; }
        }
        return result;
    }

    double IntegratedVolume(const TrimmedElementType& rElement)
    {
        double result = 0.0;
        for (const auto& rPoint : rElement.GetIntegrationPoints()) { result += rPoint.Weight() * rElement.DetJ(); }
        return result;
    }

}  // namespace

TEST_SUITE("MomentFittingComputeTestSuite")
{

    TEST_CASE("MomentFittingComputeTestSuite - Replaces Existing Integration Points")
    {
        // Verifies public Compute clears previous integration points and returns an active-domain quadrature rule.


        auto element = MakeTrimmedElement();
        constexpr IntegrationPointType sentinel({ -0.75, -0.75, -0.75 }, 123.0);
        element.GetIntegrationPoints().push_back(sentinel);

        const auto residual =
            quadrature::moment_fitting::Compute(element, { .integration_order = { 2, 2, 2 }, .residual = 1e-8 });

        REQUIRE(residual.has_value());
        CheckNear(residual.value(), 0.0, 1e-6);

        const auto& r_points = element.GetIntegrationPoints();
        CHECK_FALSE(r_points.empty());

        for (const auto& rPoint : r_points) {
            CHECK((rPoint.Weight()) > (EPS4));
            CHECK_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::parametric>(rPoint.Point()));
            CHECK_FALSE(Math::Norm(rPoint.Point() - sentinel.Point()) < 1e-12);
        }

        const double reference_volume = MeshUtilities::Volume(element.GetActiveDomainBoundaryMesh());
        CheckRelativeNear(IntegratedVolume(element), reference_volume, 1e-4);
    }

    TEST_CASE("MomentFittingComputeTestSuite - Fictitious Domain Appends Outside Points")
    {
        // Verifies fictitious-domain quadrature appends only outside-active-domain points.


        auto element_without_alpha = MakeTrimmedElement();
        auto element_with_alpha = MakeTrimmedElement();

        const auto residual_without_alpha = quadrature::moment_fitting::Compute(
            element_without_alpha, { .integration_order = { 2, 2, 2 }, .residual = 1e-8 }
        );
        const auto residual_with_alpha = quadrature::moment_fitting::Compute(
            element_with_alpha, { .integration_order = { 2, 2, 2 }, .residual = 1e-8, .fictitious_domain_alpha = 0.25 }
        );

        REQUIRE(residual_without_alpha.has_value());
        REQUIRE(residual_with_alpha.has_value());
        CheckNear(residual_without_alpha.value(), 0.0, 1e-6);
        CheckNear(residual_with_alpha.value(), 0.0, 1e-6);

        const auto& r_points_without_alpha = element_without_alpha.GetIntegrationPoints();
        const auto& r_points_with_alpha = element_with_alpha.GetIntegrationPoints();

        CHECK((CountOutsideActiveDomain(element_without_alpha)) == (0UL));
        REQUIRE((r_points_with_alpha.size()) > (r_points_without_alpha.size()));
        REQUIRE(
            (CountOutsideActiveDomain(element_with_alpha))
            == (r_points_with_alpha.size() - r_points_without_alpha.size())
        );
    }

    TEST_CASE("MomentFittingComputeTestSuite - Rejects Invalid Fictitious Domain Alpha")
    {
        // Verifies public Compute rejects fictitious-domain weight scales outside 0 < alpha <= 1.


        auto negative_alpha_element = MakeTrimmedElement();
        REQUIRE_THROWS_AS(
            quadrature::moment_fitting::Compute(
                negative_alpha_element,
                { .integration_order = { 2, 2, 2 }, .residual = 1e-8, .fictitious_domain_alpha = -0.25 }
            ),
            std::exception
        );

        auto zero_alpha_element = MakeTrimmedElement();
        REQUIRE_THROWS_AS(
            quadrature::moment_fitting::Compute(
                zero_alpha_element,
                { .integration_order = { 2, 2, 2 }, .residual = 1e-8, .fictitious_domain_alpha = 0.0 }
            ),
            std::exception
        );

        auto too_large_alpha_element = MakeTrimmedElement();
        REQUIRE_THROWS_AS(
            quadrature::moment_fitting::Compute(
                too_large_alpha_element,
                { .integration_order = { 2, 2, 2 }, .residual = 1e-8, .fictitious_domain_alpha = 1.25 }
            ),
            std::exception
        );
    }

    TEST_CASE("MomentFittingComputeTestSuite - Failed Fits Do Not Publish Weights")
    {
        const auto order_info = mf::detail::MakeIntegrationOrderInfo({ 1, 1, 1 });
        const mf::detail::IntegrationGeometry geometry{ MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }), 1.0, 0.0 };
        const mf::detail::MomentFittingProblem problem(std::vector<double>(order_info.number_of_functions, 1.0));
        for (const auto status :
             { nnls::Status::IterationLimit, nnls::Status::NumericalFailure, nnls::Status::InvalidInput }) {
            TrimmedElementType::IntegrationPointVectorType points{ IntegrationPointType({ 0.5, 0.5, 0.5 }, 7.0) };
            mf::detail::MomentFittingScratch scratch;
            const auto solve = [status](nnls::Parameters Arguments) {
                std::ranges::fill(Arguments.solution, std::numeric_limits<double>::quiet_NaN());
                return nnls::Result{ .status = status, .residual_norm = 0.0, .iterations = 1 };
            };
            if (status == nnls::Status::InvalidInput) {
                CHECK_THROWS_AS(
                    mf::detail::MomentFitting<TrimmedElementType>(
                        points, geometry, order_info, problem, scratch, solve
                    ),
                    std::exception
                );
            } else {
                const auto result = mf::detail::MomentFitting<TrimmedElementType>(
                    points, geometry, order_info, problem, scratch, solve
                );
                CHECK(result.status == status);
                CHECK(std::isinf(result.residual));
            }
            REQUIRE_FALSE(points.empty());
            CHECK_EQ(points.front().Weight(), 7.0);
        }
    }

    TEST_CASE("MomentFittingComputeTestSuite - Reduction Rollback Preserves Nonzero Support")
    {
        const auto order_info = mf::detail::MakeIntegrationOrderInfo({ 1, 1, 1 });
        const mf::detail::IntegrationGeometry geometry{ MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }), 2.0, 0.0 };
        const mf::detail::MomentFittingProblem problem(std::vector<double>(order_info.number_of_functions, 1.0));
        const auto moments = [&](const TrimmedElementType::IntegrationPointVectorType& rPoints) {
            std::vector<double> matrix, values(order_info.number_of_functions, 0.0);
            mf::detail::AssembleFittingMatrix(matrix, rPoints, geometry.bounds_param, order_info);
            for (IndexType j = 0; j < rPoints.size(); ++j) {
                for (IndexType i = 0; i < values.size(); ++i) {
                    values[i] += matrix[j * values.size() + i] * rPoints[j].Weight() * geometry.det_j;
                }
            }
            return values;
        };
        for (const auto status :
             { nnls::Status::Converged, nnls::Status::IterationLimit, nnls::Status::NumericalFailure }) {
            for (const IndexType fail_on : { IndexType{ 2 }, IndexType{ 3 } }) {
                TrimmedElementType::IntegrationPointVectorType points, expected;
                std::vector<double> expected_moments;
                for (IndexType i = 0; i < 10; ++i) { points.emplace_back(PointType{ 0.05 * i, 0.5, 0.5 }, 0.0); }
                IndexType calls = 0;
                const auto solve = [&](nnls::Parameters Arguments) {
                    if (++calls == fail_on) {
                        std::ranges::fill(
                            Arguments.solution,
                            status == nnls::Status::Converged ? 1.0 : std::numeric_limits<double>::quiet_NaN()
                        );
                        return nnls::Result{ .status = status, .residual_norm = 1.0 };
                    }
                    expected = points;
                    for (IndexType i = 0; i < points.size(); ++i) {
                        Arguments.solution[i] = static_cast<double>(i + calls);
                        if (i == 1 || i + 1 == points.size()) { Arguments.solution[i] = 0.0; }
                        if (i == 0) { Arguments.solution[i] = 1e-16; }
                        expected[i].SetWeight(Arguments.solution[i] / geometry.det_j);
                    }
                    expected_moments = moments(expected);
                    std::erase_if(expected, [](const auto& rPoint) { return rPoint.Weight() == 0.0; });
                    return nnls::Result{ .status = nnls::Status::Converged, .residual_norm = 0.0 };
                };
                const auto result = mf::detail::PointElimination<TrimmedElementType>(
                    points, geometry, order_info, problem, 1e-8, solve
                );
                CHECK(result.status == nnls::Status::Converged);
                CHECK_EQ(result.residual, 0.0);
                CHECK_EQ(calls, fail_on);
                REQUIRE_EQ(points.size(), expected.size());
                REQUIRE_FALSE(expected.empty());
                CHECK_LE(points.size(), order_info.number_of_functions);
                const auto restored_moments = moments(points);
                CHECK(restored_moments == expected_moments);
                CHECK_EQ(result.weight_stats.min_weight, expected.front().Weight());
                CHECK_EQ(result.weight_stats.min_weight_index, 0);
                CHECK_EQ(result.weight_stats.max_weight, expected.back().Weight());
                CHECK_EQ(result.weight_stats.number_of_small_weights, 1);
                for (IndexType i = 0; i < points.size(); ++i) {
                    for (IndexType axis = 0; axis < 3; ++axis) { CHECK_EQ(points[i][axis], expected[i][axis]); }
                    CHECK_EQ(points[i].Weight(), expected[i].Weight());
                }
            }
        }
    }

    TEST_CASE("MomentFittingComputeTestSuite - Failed Initial Fits Retry Then Recover Or Exclude")
    {
        for (const auto status : { nnls::Status::IterationLimit, nnls::Status::NumericalFailure }) {
            for (const bool can_recover : { false, true }) {
                auto element = MakeTrimmedElement();
                IndexType calls = 0;
                const auto solve = [&](nnls::Parameters Arguments) {
                    if (++calls == 1 || !can_recover) {
                        std::ranges::fill(Arguments.solution, std::numeric_limits<double>::quiet_NaN());
                        return nnls::Result{ .status = status, .iterations = 1 };
                    }
                    return nnls::Solve(Arguments);
                };
                const auto residual = mf::detail::ComputeWithSolver(
                    element,
                    { .integration_order = { 2, 2, 2 }, .residual = 1e-8, .fictitious_domain_alpha = 0.25 },
                    solve
                );
                CHECK_EQ(residual.has_value(), can_recover);
                if (can_recover) {
                    CHECK_GT(calls, 1);
                    CHECK(!element.GetIntegrationPoints().empty());
                    for (const auto& r_point : element.GetIntegrationPoints()) {
                        CHECK(std::isfinite(r_point.Weight()));
                        CHECK_GE(r_point.Weight(), 0.0);
                    }
                } else {
                    CHECK_EQ(calls, 3);
                    CHECK(element.GetIntegrationPoints().empty());
                }
            }
        }
    }

    TEST_CASE("MomentFittingComputeTestSuite - Failed Retries Preserve An Earlier Usable Rule")
    {
        auto element = MakeTrimmedElement();
        TrimmedElementType::IntegrationPointVectorType expected;
        IndexType calls = 0;
        const auto solve = [&](nnls::Parameters Arguments) {
            if (++calls > 1) { return nnls::Result{ .status = nnls::Status::NumericalFailure }; }
            const double rhs_norm = Math::L2Norm(Arguments.rhs);
            auto result = nnls::Solve(Arguments);
            REQUIRE(result.status == nnls::Status::Converged);
            expected = element.GetIntegrationPoints();
            for (IndexType i = 0; i < expected.size(); ++i) {
                expected[i].SetWeight(Arguments.solution[i] / element.DetJ());
            }
            std::erase_if(expected, [](const auto& rPoint) { return rPoint.Weight() == 0.0; });
            result.residual_norm = 5e-3 * rhs_norm;
            return result;
        };
        const auto residual =
            mf::detail::ComputeWithSolver(element, { .integration_order = { 2, 2, 2 }, .residual = 1e-8 }, solve);
        REQUIRE(residual.has_value());
        CHECK(std::abs(*residual - 5e-3) <= 1e-15);
        CHECK_EQ(calls, 4);
        const auto& r_points = element.GetIntegrationPoints();
        REQUIRE_EQ(r_points.size(), expected.size());
        for (IndexType i = 0; i < expected.size(); ++i) {
            for (IndexType axis = 0; axis < 3; ++axis) { CHECK_EQ(r_points[i][axis], expected[i][axis]); }
            CHECK_EQ(r_points[i].Weight(), expected[i].Weight());
        }
    }


}  // TEST_SUITE

}  // namespace queso::testing
