// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

/// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_polyhedron_3.h>

#include "math.h"

// Project includes
#include "quadrature/moment_fitting_utilities.h"
#include "modeler/modeler.h"
#include "utilities/embedding_utilities.h"
#include "quadrature/single_element.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "utilities/parameters.h"
#include "io/io_utilities.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Mesh;

namespace Testing{

BOOST_AUTO_TEST_SUITE( MomentFittingTestSuite )

BOOST_AUTO_TEST_CASE(MomentFittingP2) {
    std::cout << "Testing :: Test Moment Fitting :: Surface Integral p=2" << std::endl;

    std::array<double,3> point_a_outer = {0.0, -0.1, -0.1};
    std::array<double,3> point_b_outer = {2.1, 2.1, 3.1};

    auto geometry_outer = Modeler::make_cube_3(point_a_outer, point_b_outer);

    std::array<double,3> point_a_inner = {-0.1, 0.0, 0.0};
    std::array<double,3> point_b_inner = {2.0, 2.0, 3.0};

    auto geometry_inner = Modeler::make_cube_3(point_a_inner, point_b_inner);

    std::array<double, 3> point_A = {0.0, 0.0, 0.0};
    std::array<double, 3> point_B = {2.0, 2.0, 3.0};
    std::array<int, 3> number_of_elements = {1, 1, 1};
    std::array<int, 3> order = {2, 2, 2};

    int point_distribution_factor = 3;
    double initial_triangle_edge_length = 1;
    int minimum_number_of_triangles = 10000;
    double moment_fitting_residual = 1e-8;
    std::string integration_method = "Gauss";
    int echo_level = 0;

    Parameters param(point_A, point_B, number_of_elements, order, initial_triangle_edge_length,
        minimum_number_of_triangles, moment_fitting_residual, point_distribution_factor, integration_method, echo_level);
    param.SetUseCustomizedTrimmedPointsPositionFlag(true);

    point_B = {1.0, 1.0, 1.0};
    Element element(1, point_A, point_B, param);
    auto status = EmbeddingUtilities::ComputeIntersectionMesh( *geometry_outer, *geometry_inner, element, param);

    element.GetIntegrationPointsTrimmed().clear();

    SingleElement::Assemble(element.GetIntegrationPointsTrimmed(),
        point_A, point_B, order, param.IntegrationMethod());

    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPointsTrimmed() ){
        point.SetWeight(0.0);
    }

    SingleElement::Assemble(element.GetIntegrationPointsInside(),
        point_A, point_B, order, param.IntegrationMethod());

    // Run Moment Fitting
    MomentFitting::CreateIntegrationPointsTrimmed(element, param);

    auto& points_moment_fitting = element.GetIntegrationPointsTrimmed();
    auto& points_gauss_legendre = element.GetIntegrationPointsInside();

    double error_norm = 0.0;
    // Check if weights are similar
    for( int i = 0; i < 27; ++i){
        double weight_mf = points_moment_fitting[i].GetWeight();
        double weight_gl = points_gauss_legendre[i].GetWeight();
        double error = (weight_mf-weight_gl)/weight_gl;
        error_norm += std::pow(error,2);
        BOOST_CHECK_CLOSE_FRACTION(weight_mf, weight_gl, 1e-7);
    }
    //std::cout << "Error Norm: " << 1.0/27.0*std::sqrt(error_norm) << std::endl;

    BOOST_CHECK_LT(1.0/27.0*std::sqrt(error_norm), 1e-8);
} // End Testcase


BOOST_AUTO_TEST_CASE(MomentFittingP3) {
    std::cout << "Testing :: Test Moment Fitting :: Surface Integral p=3" << std::endl;

    std::array<double,3> point_a_outer = {0.0, -0.1, -0.1};
    std::array<double,3> point_b_outer = {1.1, 1.1, 1.1};

    auto geometry_outer = Modeler::make_cube_3(point_a_outer, point_b_outer);

    std::array<double,3> point_a_inner = {-0.1, 0.0, 0.0};
    std::array<double,3> point_b_inner = {1.0, 1.0, 1.0};

    auto geometry_inner = Modeler::make_cube_3(point_a_inner, point_b_inner);

    std::array<double, 3> point_A = {0.0, 0.0, 0.0};
    std::array<double, 3> point_B = {1.0, 1.0, 1.0};
    std::array<int, 3> number_of_elements = {1, 1, 1};
    std::array<int, 3> order = {3, 3, 3};

    int point_distribution_factor = 3;
    double initial_triangle_edge_length = 1;
    int minimum_number_of_triangles = 10000;
    double moment_fitting_residual = 1e-8;
    std::string integration_method = "Gauss";
    int echo_level = 0;

    Parameters param(point_A, point_B, number_of_elements, order, initial_triangle_edge_length,
        minimum_number_of_triangles, moment_fitting_residual, point_distribution_factor, integration_method, echo_level);
    param.SetUseCustomizedTrimmedPointsPositionFlag(true);

    Element element(1, point_A, point_B, param);
    auto status = EmbeddingUtilities::ComputeIntersectionMesh( *geometry_outer, *geometry_inner, element, param);

    element.GetIntegrationPointsTrimmed().clear();

    SingleElement::Assemble(element.GetIntegrationPointsTrimmed(),
        point_A, point_B, order, param.IntegrationMethod());

    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPointsTrimmed() ){
        point.SetWeight(0.0);
    }

    SingleElement::Assemble(element.GetIntegrationPointsInside(),
        point_A, point_B, order, param.IntegrationMethod() );

    // Run Moment Fitting
    MomentFitting::CreateIntegrationPointsTrimmed(element, param);

    auto& points_moment_fitting = element.GetIntegrationPointsTrimmed();
    auto& points_gauss_legendre = element.GetIntegrationPointsInside();

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

// BOOST_AUTO_TEST_CASE(MomentFittingVolumeIntegral) {
//     std::cout << "Testing :: Test Moment Fitting :: Volume Integral" << std::endl;

//     std::array<double,3> point_a_outer = {0.0, -0.1, -0.1};
//     std::array<double,3> point_b_outer = {1.1, 1.1, 1.1};

//     auto geometry_outer = CubeModeler::make_cube_3(point_a_outer, point_b_outer);

//     std::array<double,3> point_a_inner = {-0.1, 0.0, 0.0};
//     std::array<double,3> point_b_inner = {1.0, 1.0, 1.0};

//     auto geometry_inner = CubeModeler::make_cube_3(point_a_inner, point_b_inner);

//     std::array<double, 3> point_A = {0.0, 0.0, 0.0};
//     std::array<double, 3> point_B = {1.0, 1.0, 1.0};
//     std::array<int, 3> number_of_elements = {1, 1, 1};
//     std::array<int, 3> order = {2, 2, 2};
//     double cell_size_3D = 1;
//     double cell_size_2D = 1;
//     double cell_radius_edge_ratio = 2.0;
//     double minimum_number_of_tets = 50000;
//     double mMinimumNumberOfTriangles = 1;
//     bool reduced_points_inside_flag = false;
//     bool surface_integral_flag = false;
//     int echo_level = 0;

//     Parameters param(point_A, point_B, number_of_elements, order, cell_size_3D, cell_size_2D, cell_radius_edge_ratio,
//         minimum_number_of_tets, mMinimumNumberOfTriangles, reduced_points_inside_flag, surface_integral_flag, echo_level);

//     Element element(1, point_A, point_B);
//     auto status = EmbeddingUtilities::ComputeIntersectionMesh( *geometry_outer, *geometry_inner, element, param);

//     element.GetIntegrationPointsTrimmedReduced().clear();

//     SingleElement::Assemble(element.GetIntegrationPointsTrimmedReduced(),
//         point_A, point_B, order[0], order[1], order[2]);

//     // Make sure weights are disturbed.
//     for( auto& point : element.GetIntegrationPointsTrimmedReduced() ){
//         point.SetWeight(0.0);
//     }

//     SingleElement::Assemble(element.GetIntegrationPointsInside(),
//         point_A, point_B, order[0], order[1], order[2]);

//     // Run Moment Fitting
//     MomentFitting::ComputeReducedPointsVolumeIntegral(element, point_A, point_B);

//     auto& points_moment_fitting = element.GetIntegrationPointsTrimmedReduced();
//     auto& points_gauss_legendre = element.GetIntegrationPointsInside();

//     double error_norm = 0.0;
//     // Check if weights are similar
//     for( int i = 0; i < 27; ++i){
//         double weight_mf = points_moment_fitting[i].GetWeight();
//         double weight_gl = points_gauss_legendre[i].GetWeight();
//         double error = (weight_mf-weight_gl)/weight_gl;
//         error_norm += std::pow(error,2);
//         BOOST_CHECK_CLOSE_FRACTION(weight_mf, weight_gl, 3e-3);
//     }
//     std::cout << "Error Norm: " << error_norm << std::endl;
//     BOOST_CHECK_LT(1.0/27.0*std::sqrt(error_norm), 0.0003);
// } // End Testcase

