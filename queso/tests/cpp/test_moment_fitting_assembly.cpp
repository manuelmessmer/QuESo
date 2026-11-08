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

//// STL includes
#include <utility>

//// External includes
#include <doctest/doctest.h>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/integration_point.hpp"
#include "queso/quadrature/moment_fitting_assembly.hpp"
#include "queso/quadrature/moment_fitting_types.hpp"
#include "queso/utilities/mesh_utilities.h"
#include "queso/utilities/triangle_utilities.hpp"

// This suite tests deterministic moment-fitting assembly details: constant terms, fitting-matrix layout,
// polynomial-order dispatch and validation. Related coverage: test_moment_fitting.cpp checks NNLS solves,
// test_point_elimination.cpp checks reduced rules, and test_trimmed_domain.cpp checks STL reference data.

namespace queso::testing {

namespace {

    namespace Assembly = quadrature::moment_fitting::detail;
    using IntegrationPointType = IntegrationPoint;
    using BoundaryIntegrationPointType = BoundaryIntegrationPoint;

    void MakeAndDiscardIntegrationOrderInfo(const Vector3i& rOrder)
    { [[maybe_unused]] const auto order_info = Assembly::MakeIntegrationOrderInfo(rOrder); }

    std::vector<BoundaryIntegrationPointType> MakeBoxBoundaryIps(const BoundingBoxType& rBounds)
    {
        auto triangle_mesh = MeshUtilities::MakeMeshBox(rBounds.lower, rBounds.upper);

        std::vector<BoundaryIntegrationPointType> boundary_ips{};
        triangle_mesh.View().VisitEachTriangle<WithNormals>([&](const auto& rTriangle) {
            constexpr auto order = triangle_utilities::TriangleQuadratureOrder::Degree6;
            triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPointType>(rTriangle, order, [&](auto&& rPoint) {
                boundary_ips.push_back(std::move(rPoint));
            });
        });
        return boundary_ips;
    }

    double BoxVolume(const BoundingBoxType& rBounds)
    {
        const auto length = rBounds.upper - rBounds.lower;
        return length[0] * length[1] * length[2];
    }

    void CheckOddFirstOrderConstantTermsAreZero(const std::vector<double>& rConstantTerms)
    {
        CheckNear(rConstantTerms[1], 0.0, 1e-12);
        CheckNear(rConstantTerms[2], 0.0, 1e-12);
        CheckNear(rConstantTerms[4], 0.0, 1e-12);
    }

}  // namespace

TEST_SUITE("MomentFittingAssemblyTestSuite")
{

    TEST_CASE("MomentFittingAssemblyTestSuite - Integration Order Info Accepts Only Supported Positive Orders")
    {
        // Verifies the assembly contract: moment-fitting supports positive polynomial orders up to p=8.
        const auto order_info = Assembly::MakeIntegrationOrderInfo({ 1, 2, 3 });
        CHECK((order_info.order_u) == (1UL));
        CHECK((order_info.order_v) == (2UL));
        CHECK((order_info.order_w) == (3UL));
        CHECK((order_info.number_of_functions) == (24UL));

        const auto max_order_info = Assembly::MakeIntegrationOrderInfo({ 8, 8, 8 });
        CHECK((max_order_info.number_of_functions) == (729UL));

        REQUIRE_THROWS_AS(MakeAndDiscardIntegrationOrderInfo({ 0, 1, 1 }), std::exception);
        REQUIRE_THROWS_AS(MakeAndDiscardIntegrationOrderInfo({ 1, 9, 1 }), std::exception);
    }

    TEST_CASE("MomentFittingAssemblyTestSuite - Fitting Matrix P1 Uses Point Columns And Basis Order")
    {
        // Verifies NNLS column-major point columns and the current basis ordering u -> v -> w.
        constexpr auto bounds = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
        const std::vector<IntegrationPointType> points = {
            IntegrationPointType({ 0.25, -0.5, 0.75 }, 1.0),
            IntegrationPointType({ -0.25, 0.5, -0.75 }, 1.0),
        };

        std::vector<double> matrix{};
        Assembly::AssembleFittingMatrix(matrix, points, bounds, Assembly::MakeIntegrationOrderInfo({ 1, 1, 1 }));

        REQUIRE((matrix.size()) == (16UL));

        constexpr std::array<double, 8> first_column = { 1.0, 0.75, -0.5, -0.375, 0.25, 0.1875, -0.125, -0.09375 };
        constexpr std::array<double, 8> second_column = { 1.0, -0.75, 0.5, -0.375, -0.25, 0.1875, -0.125, 0.09375 };

        for (IndexType i = 0; i < first_column.size(); ++i) { CheckNear(matrix[i], first_column[i], 1e-12); }
        for (IndexType i = 0; i < second_column.size(); ++i) {
            CheckNear(matrix[i + first_column.size()], second_column[i], 1e-12);
        }
    }

    TEST_CASE("MomentFittingAssemblyTestSuite - Fitting Matrix Mixed Order Has Expected Size")
    {
        // Verifies that mixed directional orders use the expected tensor-product number of moment functions.
        constexpr auto bounds = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
        const std::vector<IntegrationPointType> points = {
            IntegrationPointType({ 0.0, 0.0, 0.0 }, 1.0),
            IntegrationPointType({ 0.5, -0.5, 0.25 }, 1.0),
            IntegrationPointType({ -0.5, 0.25, -0.25 }, 1.0),
        };

        std::vector<double> matrix{};
        Assembly::AssembleFittingMatrix(matrix, points, bounds, Assembly::MakeIntegrationOrderInfo({ 2, 1, 3 }));

        REQUIRE((matrix.size()) == (3UL * 2UL * 4UL * points.size()));
    }

    TEST_CASE("MomentFittingAssemblyTestSuite - Constant Terms For Symmetric Box Recover Volume And Zero Odd Moments")
    {
        // Verifies basic divergence-theorem assembly on a symmetric box without using external STL data.
        constexpr auto bounds = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
        const auto boundary_ips = MakeBoxBoundaryIps(bounds);

        const auto constant_terms =
            Assembly::ComputeConstantTerms(boundary_ips, bounds, Assembly::MakeIntegrationOrderInfo({ 1, 1, 1 }));

        REQUIRE((constant_terms.size()) == (8UL));
        CheckNear(constant_terms[0], BoxVolume(bounds), 1e-12);
        CheckOddFirstOrderConstantTermsAreZero(constant_terms);
    }

}  // TEST_SUITE

}  // namespace queso::testing
