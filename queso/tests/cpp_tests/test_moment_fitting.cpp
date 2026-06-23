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
#include "math.h"
// Project includes
#include "queso/includes/checks.hpp"
#include "queso/containers/untrimmed_element.hpp"
#include "queso/quadrature/single_element.hpp"
#include "queso/tests/cpp_tests/class_testers/trimmed_element_tester.hpp"
#include "queso/utilities/mesh_utilities.h"
#include "queso/utilities/triangle_utilities.hpp"
#include "queso/containers/boundary_integration_point.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( MomentFittingTestSuite )

BOOST_AUTO_TEST_CASE(MomentFittingP2) {
    QuESo_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=2" << std::endl;

    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef UntrimmedElement<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    const auto bounds_xyz = MakeBox({0, 0, 0}, {1, 1, 3.0});
    const auto bounds_uvw = MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
    ElementType element(1, ElementBounds{bounds_xyz, bounds_uvw});

    // Construct cube over domian.
    PointType point_a_domain = {0.0, 0.0, 0.0};
    PointType point_b_domain = {1.0, 1.0, 3.0};
    auto triangle_mesh = MeshUtilities::MakeMeshBox(point_a_domain, point_b_domain);

    ElementType::BoundaryIntegrationPointVectorType boundary_ips{};
    triangle_mesh.View().VisitEachTriangle<WithNormals>([&](const auto& triangle) {
        const IndexType method = 3; // This will create 6 points per triangle.
        auto new_points = TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPointType>(triangle, method);
        boundary_ips.insert(boundary_ips.end(), new_points.begin(), new_points.end());
    });

    const Vector3i polynomial_order = {2, 2, 2};
    const IntegrationMethod integration_method = IntegrationMethod::gauss;

    // Distribtue Gauss points within element.
    element.GetIntegrationPoints().clear();
    QuadratureSingleElement<ElementType>::AssembleIPs(element, polynomial_order, integration_method);
    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    // Run Moment Fitting
    std::vector<double> constant_terms{};
    QuadratureTrimmedElementTester<ElementType>::ComputeConstantTerms(constant_terms, boundary_ips, element, polynomial_order);
    QuadratureTrimmedElementTester<ElementType>::MomentFitting(constant_terms, element.GetIntegrationPoints(), element, polynomial_order);
    auto& points_moment_fitting = element.GetIntegrationPoints();

    // Get Gauss points as reference
    ElementType::IntegrationPointVectorType points_gauss_legendre{};
    QuadratureSingleElement<ElementType>::AssembleIPs(points_gauss_legendre, element.GetCellBounds<CoordinateSpace::parametric>().lower,
        element.GetCellBounds<CoordinateSpace::parametric>().upper, {2, 2, 2});

    double error_norm = 0.0;
    // Check if weights are similar
    for( int i = 0; i < 27; ++i){
        double weight_mf = points_moment_fitting[i].Weight();
        double weight_gl = points_gauss_legendre[i].Weight();
        double error = (weight_mf-weight_gl)/weight_gl;
        error_norm += std::pow(error,2);
        QuESo_CHECK_RELATIVE_NEAR(weight_mf, weight_gl, 1e-12);
    }

    QuESo_CHECK_LT(1.0/27.0*std::sqrt(error_norm), 1e-10);
} // End Testcase


BOOST_AUTO_TEST_CASE(MomentFittingP3) {
    QuESo_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=3" << std::endl;

    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef UntrimmedElement<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    const auto bounds_xyz = MakeBox({0, 0, 0}, {2.0, 2.0, 1.0});
    const auto bounds_uvw = MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
    ElementType element(1, ElementBounds{bounds_xyz, bounds_uvw});

    // Construct cube over domian.
    PointType point_a_domain = {0.0, 0.0, 0.0};
    PointType point_b_domain = {2.0, 2.0, 1.0};
    auto triangle_mesh = MeshUtilities::MakeMeshBox(point_a_domain, point_b_domain);

    MeshUtilities::Refine(triangle_mesh, 500 );
    ElementType::BoundaryIntegrationPointVectorType boundary_ips{};
    triangle_mesh.View().VisitEachTriangle<WithNormals>([&](const auto& triangle) {
        const IndexType method = 3; // This will create 6 points per triangle.
        auto new_points = TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPointType>(triangle, method);
        boundary_ips.insert(boundary_ips.end(), new_points.begin(), new_points.end());
    });

    const Vector3i polynomial_order = {3, 3, 3};
    const IntegrationMethod integration_method = IntegrationMethod::gauss;

    // Distribtue Gauss points within element.
    element.GetIntegrationPoints().clear();
    QuadratureSingleElement<ElementType>::AssembleIPs(element, polynomial_order, integration_method);
    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    // Run Moment Fitting
    std::vector<double> constant_terms{};
    QuadratureTrimmedElementTester<ElementType>::ComputeConstantTerms(constant_terms, boundary_ips, element, polynomial_order);
    QuadratureTrimmedElementTester<ElementType>::MomentFitting(constant_terms, element.GetIntegrationPoints(), element, polynomial_order);
    auto& points_moment_fitting = element.GetIntegrationPoints();

    // Get Gauss points as reference
    ElementType::IntegrationPointVectorType points_gauss_legendre{};
    QuadratureSingleElement<ElementType>::AssembleIPs(points_gauss_legendre, element.GetCellBounds<CoordinateSpace::parametric>().lower,
        element.GetCellBounds<CoordinateSpace::parametric>().upper, {3, 3, 3});

    double error_norm = 0.0;
    // Check if weights are similar
    for( int i = 0; i < 64; ++i){
        double weight_mf = points_moment_fitting[i].Weight();
        double weight_gl = points_gauss_legendre[i].Weight();
        double error = (weight_mf-weight_gl)/weight_gl;
        error_norm += std::pow(error,2);
        QuESo_CHECK_RELATIVE_NEAR(weight_mf, weight_gl, 1e-6);
    }
    //QuESo_INFO << "Error Norm: " << 1.0/64.0*std::sqrt(error_norm) << std::endl;

    QuESo_CHECK_LT(1.0/64.0*std::sqrt(error_norm), 1e-7);
} // End Testcase

BOOST_AUTO_TEST_CASE(MomentFittingP4) {
    QuESo_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=4" << std::endl;

    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef UntrimmedElement<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    const auto bounds_xyz = MakeBox({0, 0, 0}, {2.0, 2.0, 1.0});
    const auto bounds_uvw = MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
    ElementType element(1, ElementBounds{bounds_xyz, bounds_uvw});

    // Construct cube over domian.
    PointType point_a_domain = {0.0, 0.0, 0.0};
    PointType point_b_domain = {2.0, 2.0, 1.0};
    auto triangle_mesh = MeshUtilities::MakeMeshBox(point_a_domain, point_b_domain);

    MeshUtilities::Refine(triangle_mesh, 2000 );
    ElementType::BoundaryIntegrationPointVectorType boundary_ips{};
    triangle_mesh.View().VisitEachTriangle<WithNormals>([&](const auto& triangle) {
        const IndexType method = 3; // This will create 6 points per triangle.
        auto new_points = TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPointType>(triangle, method);
        boundary_ips.insert(boundary_ips.end(), new_points.begin(), new_points.end());
    });

    const Vector3i polynomial_order = {4, 4, 4};
    const IntegrationMethod integration_method = IntegrationMethod::gauss;

    // Distribtue Gauss points within element.
    element.GetIntegrationPoints().clear();
    QuadratureSingleElement<ElementType>::AssembleIPs(element, polynomial_order, integration_method);
    // Make sure weights are disturbed.param
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    // Run Moment Fitting
    std::vector<double> constant_terms{};
    QuadratureTrimmedElementTester<ElementType>::ComputeConstantTerms(constant_terms, boundary_ips, element, polynomial_order);
    QuadratureTrimmedElementTester<ElementType>::MomentFitting(constant_terms, element.GetIntegrationPoints(), element, polynomial_order);
    auto& points_moment_fitting = element.GetIntegrationPoints();

    // Get Gauss points as reference
    ElementType::IntegrationPointVectorType points_gauss_legendre{};
    QuadratureSingleElement<ElementType>::AssembleIPs(points_gauss_legendre, element.GetCellBounds<CoordinateSpace::parametric>().lower,
        element.GetCellBounds<CoordinateSpace::parametric>().upper, polynomial_order);

    double error_norm = 0.0;
    // Check if weights are similar
    for( int i = 0; i < 125; ++i){
        double weight_mf = points_moment_fitting[i].Weight();
        double weight_gl = points_gauss_legendre[i].Weight();
        double error = (weight_mf-weight_gl)/weight_gl;
        error_norm += std::pow(error,2);
        QuESo_CHECK_RELATIVE_NEAR(weight_mf, weight_gl, 1e-6);
    }
    //QuESo_INFO << "Error Norm: " << 1.0/64.0*std::sqrt(error_norm) << std::endl;

    QuESo_CHECK_LT(1.0/125.0*std::sqrt(error_norm), 1e-7);
} // End Testcase

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso
