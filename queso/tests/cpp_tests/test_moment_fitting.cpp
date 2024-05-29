// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// STL includes
#include "math.h"
// Project includes
#include "includes/checks.hpp"
#include "embedding/brep_operator.h"
#include "quadrature/trimmed_element.hpp"
#include "quadrature/single_element.hpp"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "includes/parameters.h"
#include "utilities/mapping_utilities.h"
#include "utilities/mesh_utilities.h"
#include "io/io_utilities.h"
#include "tests/cpp_tests/class_testers/trimmed_element_tester.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( MomentFittingTestSuite )

BOOST_AUTO_TEST_CASE(MomentFittingP2) {
    QuESo_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=2" << std::endl;

    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    Parameters parameters( {Component("lower_bound_xyz", PointType{0.0, 0.0, 0.0}),
                            Component("upper_bound_xyz", PointType{2.0, 2.0, 3.0}),
                            Component("number_of_elements", Vector3i{1, 1, 1}),
                            Component("polynomial_order", Vector3i{2, 2, 2}),
                            Component("moment_fitting_residual", 1e-8),
                            Component("min_num_boundary_triangles", 5UL),
                            Component("integration_method", IntegrationMethod::Gauss),
                            Component("use_customized_trimmed_points", false) });

    ElementType element(1, MakeBox({0, 0, 0}, {1, 1, 3.0}), MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}));

    // Construct cube over domian.
    PointType point_a_domain = {0.0, 0.0, 0.0};
    PointType point_b_domain = {1.0, 1.0, 3.0};
    auto p_triangle_mesh = MeshUtilities::pGetCuboid(point_a_domain, point_b_domain);

    auto p_boundary_ips = MakeUnique<ElementType::BoundaryIntegrationPointVectorType>();
    for( IndexType triangle_id = 0; triangle_id < p_triangle_mesh->NumOfTriangles(); ++triangle_id ) {
            IndexType method = 3; // This will create 6 points per triangle.
            auto p_new_points = p_triangle_mesh->pGetIPsGlobal<BoundaryIntegrationPointType>(triangle_id, method);
            p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
    }

    const Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    const IntegrationMethod integration_method = parameters.IntegrationMethod();

    // Distribtue Gauss points within element.
    element.GetIntegrationPoints().clear();
    QuadratureSingleElement<ElementType>::AssembleIPs(element, polynomial_order, integration_method);
    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    // Run Moment Fitting
    std::vector<double> constant_terms{};
    QuadratureTrimmedElementTester<ElementType>::ComputeConstantTerms(constant_terms, p_boundary_ips, element, polynomial_order);
    QuadratureTrimmedElementTester<ElementType>::MomentFitting(constant_terms, element.GetIntegrationPoints(), element, polynomial_order);
    auto& points_moment_fitting = element.GetIntegrationPoints();

    // Get Gauss points as reference
    ElementType::IntegrationPointVectorType points_gauss_legendre{};
    QuadratureSingleElement<ElementType>::AssembleIPs(points_gauss_legendre, element.GetBoundsUVW().first,
        element.GetBoundsUVW().second, {2, 2, 2});

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
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    Parameters parameters( {Component("lower_bound_xyz", PointType{0.0, 0.0, 0.0}),
                            Component("upper_bound_xyz", PointType{2.0, 2.0, 1.0}),
                            Component("number_of_elements", Vector3i{1, 1, 1}),
                            Component("polynomial_order", Vector3i{3, 3, 3}),
                            Component("moment_fitting_residual", 1e-8),
                            Component("min_num_boundary_triangles", 500UL),
                            Component("integration_method", IntegrationMethod::Gauss),
                            Component("use_customized_trimmed_points", false) });

    ElementType element(1, MakeBox({0, 0, 0}, {2.0, 2.0, 1.0}), MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}));

    // Construct cube over domian.
    PointType point_a_domain = {0.0, 0.0, 0.0};
    PointType point_b_domain = {2.0, 2.0, 1.0};
    auto p_triangle_mesh = MeshUtilities::pGetCuboid(point_a_domain, point_b_domain);

    MeshUtilities::Refine(*p_triangle_mesh, parameters.MinimumNumberOfTriangles() );
    auto p_boundary_ips = MakeUnique<ElementType::BoundaryIntegrationPointVectorType>();
    for( IndexType triangle_id = 0; triangle_id < p_triangle_mesh->NumOfTriangles(); ++triangle_id ) {
            IndexType method = 3; // This will create 6 points per triangle.
            auto p_new_points = p_triangle_mesh->pGetIPsGlobal<BoundaryIntegrationPoint>(triangle_id, method);
            p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
    }

    const Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    const IntegrationMethod integration_method = parameters.IntegrationMethod();

    // Distribtue Gauss points within element.
    element.GetIntegrationPoints().clear();
    QuadratureSingleElement<ElementType>::AssembleIPs(element, polynomial_order, integration_method);
    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    // Run Moment Fitting
    std::vector<double> constant_terms{};
    QuadratureTrimmedElementTester<ElementType>::ComputeConstantTerms(constant_terms, p_boundary_ips, element, polynomial_order);
    QuadratureTrimmedElementTester<ElementType>::MomentFitting(constant_terms, element.GetIntegrationPoints(), element, polynomial_order);
    auto& points_moment_fitting = element.GetIntegrationPoints();

    // Get Gauss points as reference
    ElementType::IntegrationPointVectorType points_gauss_legendre{};
    QuadratureSingleElement<ElementType>::AssembleIPs(points_gauss_legendre, element.GetBoundsUVW().first,
        element.GetBoundsUVW().second, {3, 3, 3});

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
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    Parameters parameters( {Component("lower_bound_xyz", PointType{0.0, 0.0, 0.0}),
                            Component("upper_bound_xyz", PointType{2.0, 2.0, 1.0}),
                            Component("number_of_elements", Vector3i{1, 1, 1}),
                            Component("polynomial_order", Vector3i{4, 4, 4}),
                            Component("min_num_boundary_triangles", 2000UL),
                            Component("integration_method", IntegrationMethod::Gauss),
                            Component("use_customized_trimmed_points", false) });

     ElementType element(1, MakeBox({0, 0, 0}, {2.0, 2.0, 1.0}), MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}));

    // Construct cube over domian.
    PointType point_a_domain = {0.0, 0.0, 0.0};
    PointType point_b_domain = {2.0, 2.0, 1.0};
    auto p_triangle_mesh = MeshUtilities::pGetCuboid(point_a_domain, point_b_domain);

    MeshUtilities::Refine(*p_triangle_mesh, parameters.MinimumNumberOfTriangles() );
    auto p_boundary_ips = MakeUnique<ElementType::BoundaryIntegrationPointVectorType>();
    for( IndexType triangle_id = 0; triangle_id < p_triangle_mesh->NumOfTriangles(); ++triangle_id ) {
            IndexType method = 3; // This will create 6 points per triangle.
            auto p_new_points = p_triangle_mesh->pGetIPsGlobal<BoundaryIntegrationPointType>(triangle_id, method);
            p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
    }

    const Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    const IntegrationMethod integration_method = parameters.IntegrationMethod();

    // Distribtue Gauss points within element.
    element.GetIntegrationPoints().clear();
    QuadratureSingleElement<ElementType>::AssembleIPs(element, polynomial_order, integration_method);
    // Make sure weights are disturbed.param
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    // Run Moment Fitting
    std::vector<double> constant_terms{};
    QuadratureTrimmedElementTester<ElementType>::ComputeConstantTerms(constant_terms, p_boundary_ips, element, polynomial_order);
    QuadratureTrimmedElementTester<ElementType>::MomentFitting(constant_terms, element.GetIntegrationPoints(), element, polynomial_order);
    auto& points_moment_fitting = element.GetIntegrationPoints();

    // Get Gauss points as reference
    ElementType::IntegrationPointVectorType points_gauss_legendre{};
    QuadratureSingleElement<ElementType>::AssembleIPs(points_gauss_legendre, element.GetBoundsUVW().first,
        element.GetBoundsUVW().second, polynomial_order);

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
