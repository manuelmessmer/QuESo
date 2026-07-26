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
#include <boost/test/unit_test.hpp>

//// Project includes
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/integration_point.hpp"
#include "queso/includes/checks.hpp"
#include "queso/quadrature/moment_fitting_assembly.hpp"
#include "queso/quadrature/moment_fitting_types.hpp"
#include "queso/utilities/mesh_utilities.h"
#include "queso/utilities/triangle_utilities.hpp"

// This suite tests deterministic moment-fitting assembly details: constant terms, fitting-matrix layout,
// polynomial-order dispatch and validation. Related coverage: test_moment_fitting.cpp checks NNLS solves,
// test_point_elimination.cpp checks reduced rules, and test_trimmed_domain.cpp checks STL reference data.

namespace queso::Testing {

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
            constexpr IndexType method = 3;
            auto new_points = TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPointType>(rTriangle, method);
            boundary_ips.insert(boundary_ips.end(), new_points.begin(), new_points.end());
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
        QuESo_CHECK_NEAR(rConstantTerms[1], 0.0, 1e-12);
        QuESo_CHECK_NEAR(rConstantTerms[2], 0.0, 1e-12);
        QuESo_CHECK_NEAR(rConstantTerms[4], 0.0, 1e-12);
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(MomentFittingAssemblyTestSuite)

BOOST_AUTO_TEST_CASE(IntegrationOrderInfoAcceptsOnlySupportedPositiveOrders)
{
    // Verifies the assembly contract: moment-fitting supports positive polynomial orders up to p=8.
    const auto order_info = Assembly::MakeIntegrationOrderInfo({ 1, 2, 3 });
    QuESo_CHECK_EQUAL(order_info.order_u, 1UL);
    QuESo_CHECK_EQUAL(order_info.order_v, 2UL);
    QuESo_CHECK_EQUAL(order_info.order_w, 3UL);
    QuESo_CHECK_EQUAL(order_info.number_of_functions, 24UL);

    const auto max_order_info = Assembly::MakeIntegrationOrderInfo({ 8, 8, 8 });
    QuESo_CHECK_EQUAL(max_order_info.number_of_functions, 729UL);

    BOOST_REQUIRE_THROW(MakeAndDiscardIntegrationOrderInfo({ 0, 1, 1 }), std::exception);
    BOOST_REQUIRE_THROW(MakeAndDiscardIntegrationOrderInfo({ 1, 9, 1 }), std::exception);
}

BOOST_AUTO_TEST_CASE(FittingMatrixP1UsesPointColumnsAndBasisOrder)
{
    // Verifies NNLS column-major point columns and the current basis ordering u -> v -> w.
    constexpr auto bounds = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
    const std::vector<IntegrationPointType> points = {
        IntegrationPointType({ 0.25, -0.5, 0.75 }, 1.0),
        IntegrationPointType({ -0.25, 0.5, -0.75 }, 1.0),
    };

    NNLS::MatrixType matrix{};
    Assembly::AssembleFittingMatrix(matrix, points, bounds, Assembly::MakeIntegrationOrderInfo({ 1, 1, 1 }));

    QuESo_CHECK_EQUAL(matrix.size(), 16UL);

    constexpr std::array<double, 8> first_column = { 1.0, 0.75, -0.5, -0.375, 0.25, 0.1875, -0.125, -0.09375 };
    constexpr std::array<double, 8> second_column = { 1.0, -0.75, 0.5, -0.375, -0.25, 0.1875, -0.125, 0.09375 };

    for (IndexType i = 0; i < first_column.size(); ++i) { QuESo_CHECK_NEAR(matrix[i], first_column[i], 1e-12); }
    for (IndexType i = 0; i < second_column.size(); ++i) {
        QuESo_CHECK_NEAR(matrix[i + first_column.size()], second_column[i], 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(FittingMatrixMixedOrderHasExpectedSize)
{
    // Verifies that mixed directional orders use the expected tensor-product number of moment functions.
    constexpr auto bounds = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
    const std::vector<IntegrationPointType> points = {
        IntegrationPointType({ 0.0, 0.0, 0.0 }, 1.0),
        IntegrationPointType({ 0.5, -0.5, 0.25 }, 1.0),
        IntegrationPointType({ -0.5, 0.25, -0.25 }, 1.0),
    };

    NNLS::MatrixType matrix{};
    Assembly::AssembleFittingMatrix(matrix, points, bounds, Assembly::MakeIntegrationOrderInfo({ 2, 1, 3 }));

    QuESo_CHECK_EQUAL(matrix.size(), 3UL * 2UL * 4UL * points.size());
}

BOOST_AUTO_TEST_CASE(ConstantTermsForSymmetricBoxRecoverVolumeAndZeroOddMoments)
{
    // Verifies basic divergence-theorem assembly on a symmetric box without using external STL data.
    constexpr auto bounds = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
    const auto boundary_ips = MakeBoxBoundaryIps(bounds);

    const auto constant_terms =
        Assembly::ComputeConstantTerms(boundary_ips, bounds, Assembly::MakeIntegrationOrderInfo({ 1, 1, 1 }));

    QuESo_CHECK_EQUAL(constant_terms.size(), 8UL);
    QuESo_CHECK_NEAR(constant_terms[0], BoxVolume(bounds), 1e-12);
    CheckOddFirstOrderConstantTermsAreZero(constant_terms);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
