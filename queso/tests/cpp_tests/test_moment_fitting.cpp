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

//// STL includes
#include <cmath>

//// Project includes
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/untrimmed_element.hpp"
#include "queso/includes/checks.hpp"
#include "queso/quadrature/moment_fitting.hpp"
#include "queso/quadrature/tensor_product.hpp"
#include "queso/utilities/mesh_utilities.h"
#include "queso/utilities/triangle_utilities.hpp"

// This suite tests the moment-fitting NNLS solve on simple untrimmed boxes by verifying that fitted weights reproduce
// tensor-product Gauss weights. Related coverage: test_moment_fitting_assembly.cpp checks matrix/RHS assembly details,
// and test_point_elimination.cpp checks reduced quadrature rules on trimmed geometries.

namespace queso::Testing {

namespace {

    using IntegrationPointType = IntegrationPoint;
    using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
    using ElementType = UntrimmedElement<IntegrationPointType, BoundaryIntegrationPointType>;

    ElementType::BoundaryIntegrationPointVectorType
        MakeBoxBoundaryIps(const BoundingBoxType& rBounds, IndexType MinNumberOfBoundaryTriangles)
    {
        auto triangle_mesh = MeshUtilities::MakeMeshBox(rBounds.lower, rBounds.upper);
        if (MinNumberOfBoundaryTriangles > triangle_mesh.NumOfTriangles()) {
            MeshUtilities::Refine(
                triangle_mesh,
                MinNumberOfBoundaryTriangles,
                GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 }).ZeroArea()
            );
        }

        ElementType::BoundaryIntegrationPointVectorType boundary_ips{};
        triangle_mesh.View().VisitEachTriangle<WithNormals>([&](const auto& rTriangle) {
            constexpr IndexType method = 3;
            auto new_points = TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPointType>(rTriangle, method);
            boundary_ips.insert(boundary_ips.end(), new_points.begin(), new_points.end());
        });
        return boundary_ips;
    }

    void ResetWeights(ElementType::IntegrationPointVectorType& rPoints)
    {
        for (auto& rPoint : rPoints) { rPoint.SetWeight(0.0); }
    }

    double RelativeWeightErrorNorm(
        const ElementType::IntegrationPointVectorType& rMomentFittingPoints,
        const ElementType::IntegrationPointVectorType& rReferencePoints
    )
    {
        double error_norm = 0.0;
        for (IndexType i = 0; i < rReferencePoints.size(); ++i) {
            const double weight_mf = rMomentFittingPoints[i].Weight();
            const double weight_ref = rReferencePoints[i].Weight();
            const double error = (weight_mf - weight_ref) / weight_ref;
            error_norm += error * error;
        }
        return std::sqrt(error_norm / rReferencePoints.size());
    }

    void CheckMomentFittingRecoversTensorProductWeights(
        const BoundingBoxType& rBoundsXYZ,
        const BoundingBoxType& rBoundsUVW,
        const Vector3i& rPolynomialOrder,
        IndexType MinNumberOfBoundaryTriangles,
        double WeightTolerance,
        double ErrorNormTolerance
    )
    {
        ElementType element(1, ElementBounds{ rBoundsXYZ, rBoundsUVW });
        const auto boundary_ips = MakeBoxBoundaryIps(rBoundsXYZ, MinNumberOfBoundaryTriangles);

        quadrature::tensor_product::Compute(
            element, { .integration_order = rPolynomialOrder, .method = IntegrationMethod::gauss }
        );
        ResetWeights(element.GetIntegrationPoints());

        quadrature::moment_fitting::detail::MomentFittingScratch scratch{};
        auto constant_terms = quadrature::moment_fitting::detail::ComputeConstantTerms(
            boundary_ips,
            element.GetCellBounds<CoordinateSpace::global>(),
            quadrature::moment_fitting::detail::MakeIntegrationOrderInfo(rPolynomialOrder)
        );
        const quadrature::moment_fitting::detail::MomentFittingProblem problem(std::move(constant_terms));
        const quadrature::moment_fitting::detail::IntegrationGeometry geometry{
            element.GetCellBounds<CoordinateSpace::parametric>(), element.DetJ()
        };
        quadrature::moment_fitting::detail::MomentFitting<ElementType>(
            element.GetIntegrationPoints(),
            geometry,
            quadrature::moment_fitting::detail::MakeIntegrationOrderInfo(rPolynomialOrder),
            problem,
            scratch
        );
        const auto& r_points_moment_fitting = element.GetIntegrationPoints();

        ElementType::IntegrationPointVectorType reference_points{};
        quadrature::tensor_product::Compute(
            reference_points,
            element.GetCellBounds<CoordinateSpace::parametric>(),
            { .integration_order = rPolynomialOrder }
        );

        QuESo_CHECK_EQUAL(r_points_moment_fitting.size(), reference_points.size());
        for (IndexType i = 0; i < reference_points.size(); ++i) {
            QuESo_CHECK_RELATIVE_NEAR(
                r_points_moment_fitting[i].Weight(), reference_points[i].Weight(), WeightTolerance
            );
        }
        QuESo_CHECK_LT(RelativeWeightErrorNorm(r_points_moment_fitting, reference_points), ErrorNormTolerance);
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(MomentFittingTestSuite)

BOOST_AUTO_TEST_CASE(MomentFittingP1)
{
    // Verifies that a linear moment-fitting solve reproduces the tensor-product Gauss rule on a box.
    QuESo_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=1" << std::endl;
    constexpr auto bounds_xyz = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 3.0 });
    constexpr auto bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
    constexpr Vector3i polynomial_order{ 1, 1, 1 };
    CheckMomentFittingRecoversTensorProductWeights(bounds_xyz, bounds_uvw, polynomial_order, 0, 1e-12, 1e-10);
}

BOOST_AUTO_TEST_CASE(MomentFittingP2)
{
    // Verifies that a quadratic moment-fitting solve reproduces the tensor-product Gauss rule on a box.
    QuESo_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=2" << std::endl;
    constexpr auto bounds_xyz = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 3.0 });
    constexpr auto bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
    constexpr Vector3i polynomial_order{ 2, 2, 2 };
    CheckMomentFittingRecoversTensorProductWeights(bounds_xyz, bounds_uvw, polynomial_order, 0, 1e-12, 1e-10);
}

BOOST_AUTO_TEST_CASE(MomentFittingP3)
{
    // Verifies that a cubic moment-fitting solve reproduces the tensor-product Gauss rule on a box.
    QuESo_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=3" << std::endl;
    constexpr auto bounds_xyz = MakeBox({ 0.0, 0.0, 0.0 }, { 2.0, 2.0, 1.0 });
    constexpr auto bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
    constexpr Vector3i polynomial_order{ 3, 3, 3 };
    CheckMomentFittingRecoversTensorProductWeights(bounds_xyz, bounds_uvw, polynomial_order, 500, 1e-6, 1e-7);
}

BOOST_AUTO_TEST_CASE(MomentFittingP4)
{
    // Verifies that a quartic moment-fitting solve reproduces the tensor-product Gauss rule on a box.
    QuESo_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=4" << std::endl;
    constexpr auto bounds_xyz = MakeBox({ 0.0, 0.0, 0.0 }, { 2.0, 2.0, 1.0 });
    constexpr auto bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
    constexpr Vector3i polynomial_order{ 4, 4, 4 };
    CheckMomentFittingRecoversTensorProductWeights(bounds_xyz, bounds_uvw, polynomial_order, 2000, 1e-6, 1e-7);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
