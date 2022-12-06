// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// STL includes
#include "math.h"
// Project includes
#include "embedding/brep_operator_factory.h"
#include "quadrature/moment_fitting_utilities.h"
#include "quadrature/single_element.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "utilities/parameters.h"
#include "utilities/mapping_utilities.h"
#include "io/io_utilities.h"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( MomentFittingTestSuite )

BOOST_AUTO_TEST_CASE(MomentFittingP2) {
    std::cout << "Testing :: Test Moment Fitting :: Surface Integral p=2" << std::endl;


    Parameters parameters( {Component("lower_bound", PointType(0.0, 0.0, 0.0)),
                            Component("upper_bound", PointType(2.0, 2.0, 3.0)),
                            Component("number_of_knot_spans", Vector3i(1, 1, 1)),
                            Component("polynomial_order", Vector3i(2, 2, 2)),
                            Component("moment_fitting_residual", 1e-8),
                            Component("min_num_boundary_triangles", 20000UL),
                            Component("init_point_distribution_factor", 3UL),
                            Component("use_customized_trimmed_points", true) });

    Element element(1, {0, 0, 0}, {0.5, 0.5, 1}, parameters);

    PointType point_a_domain = {0.0, -0.1, -0.1};
    PointType point_b_domain = {2.1, 2.1, 3.1};
    auto p_triangle_mesh = TriangleMesh::MakeCuboid(point_a_domain, point_b_domain);

    element.GetIntegrationPoints().clear();
    SingleElement::AssembleIPs(element, parameters);

    // Set positions of trimmed points at Gauss Legendre Points.
    element.GetIntegrationPoints() = element.GetIntegrationPoints();
    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    // Create BRepOperator
    auto p_brep_operator = BRepOperatorFactory::New(*p_triangle_mesh);

    // Get Trimmed Domain
    auto p_trimmed_domain = p_brep_operator->GetTrimmedDomain(
        element.GetLowerBound(), element.GetUpperBound(), parameters);

    // Check if ptr is not nullptr
    BOOST_CHECK( p_trimmed_domain );
    element.SetIsTrimmed(true);
    element.pSetTrimmedDomain(p_trimmed_domain);

    // Run Moment Fitting
    MomentFitting::CreateIntegrationPointsTrimmed(element, parameters);

    auto& points_moment_fitting = element.GetIntegrationPoints();
    auto& points_gauss_legendre = element.GetIntegrationPoints();

    double error_norm = 0.0;
    // Check if weights are similar
    for( int i = 0; i < 27; ++i){
        double weight_mf = points_moment_fitting[i].GetWeight();
        double weight_gl = points_gauss_legendre[i].GetWeight();
        double error = (weight_mf-weight_gl)/weight_gl;
        error_norm += std::pow(error,2);
        BOOST_CHECK_CLOSE_FRACTION(weight_mf, weight_gl, 1e-5);
    }
    //std::cout << "Error Norm: " << 1.0/27.0*std::sqrt(error_norm) << std::endl;

    BOOST_CHECK_LT(1.0/27.0*std::sqrt(error_norm), 1e-7);
} // End Testcase


BOOST_AUTO_TEST_CASE(MomentFittingP3) {
    std::cout << "Testing :: Test Moment Fitting :: Surface Integral p=3" << std::endl;

    Parameters parameters( {Component("lower_bound", PointType(0.0, 0.0, 0.0)),
                            Component("upper_bound", PointType(2.0, 2.0, 1.0)),
                            Component("number_of_knot_spans", Vector3i(1, 1, 1)),
                            Component("polynomial_order", Vector3i(3, 3, 3)),
                            Component("moment_fitting_residual", 1e-8),
                            Component("min_num_boundary_triangles", 10000UL),
                            Component("init_point_distribution_factor", 3UL),
                            Component("use_customized_trimmed_points", true) });

    Element element(1, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, parameters);

    PointType point_a_domain = {0.0, -0.1, -0.1};
    PointType point_b_domain = {2.1, 2.1, 1.1};
    auto p_triangle_mesh = TriangleMesh::MakeCuboid(point_a_domain, point_b_domain);

    element.GetIntegrationPoints().clear();
    SingleElement::AssembleIPs(element, parameters);

    // Set positions of trimmed points at Gauss Legendre Points.
    element.GetIntegrationPoints() = element.GetIntegrationPoints();
    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    // Create BRepOperator
    auto p_brep_operator = BRepOperatorFactory::New(*p_triangle_mesh);

    // Get Trimmed Domain
    auto p_trimmed_domain = p_brep_operator->GetTrimmedDomain(
        element.GetLowerBound(), element.GetUpperBound(), parameters);
    // Check if ptr is not nullptr
    BOOST_CHECK( p_trimmed_domain );
    element.SetIsTrimmed(true);
    element.pSetTrimmedDomain(p_trimmed_domain);

    // Run Moment Fitting
    MomentFitting::CreateIntegrationPointsTrimmed(element, parameters);

    auto& points_moment_fitting = element.GetIntegrationPoints();
    auto& points_gauss_legendre = element.GetIntegrationPoints();

    double error_norm = 0.0;
    // Check if weights are similar
    for( int i = 0; i < 64; ++i){
        double weight_mf = points_moment_fitting[i].GetWeight();
        double weight_gl = points_gauss_legendre[i].GetWeight();
        double error = (weight_mf-weight_gl)/weight_gl;
        error_norm += std::pow(error,2);
        BOOST_CHECK_CLOSE_FRACTION(weight_mf, weight_gl, 1e-6);
    }
    //std::cout << "Error Norm: " << 1.0/64.0*std::sqrt(error_norm) << std::endl;

    BOOST_CHECK_LT(1.0/64.0*std::sqrt(error_norm), 1e-7);
} // End Testcase

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra
