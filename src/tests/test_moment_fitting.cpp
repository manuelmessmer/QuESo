// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "BaseClassModule"

// External includes
#include <boost/test/unit_test.hpp>
#include "math.h"

// Project includes
#include "utilities/inside_test.h"
#include "utilities/moment_fitting_utilities.h"
#include "modeler/cube_modeler.h"
#include "utilities/mesh_utilities.h"
#include "utilities/integration_point_utilities.h"
#include "utilities/integration_points/integration_points_factory.h"
#include "utilities/parameters.h"

namespace Testing{

BOOST_AUTO_TEST_CASE(MomentFittingTestLegendrePolynomials1) {
    std::cout << "Testing :: Test Moment Fitting :: Legendre Polynomials1" << std::endl;
    for(int order = 1; order <= 9; ++order){
        for( int order2 = 1; order2 <= 9; ++order2){
            if( order != order2){
                auto ips_1 = IntegrationPointFactory::GetIntegrationPoints(order, IntegrationPointFactory::IntegrationMethod::Gauss);
                auto ips_2 = IntegrationPointFactory::GetIntegrationPoints(order2, IntegrationPointFactory::IntegrationMethod::Gauss);
                double numerical_integral = 0.0;
                for( auto& point1 : ips_1){
                    for( auto& point2 : ips_2){
                        double position1 = point1[0]* 0.2 + 0.1;
                        double position2 = point2[0]* 0.2 + 0.1;
                        numerical_integral += MomentFitting::f_x(position1, order-1, 0.1, 0.3)*point1[1]
                            *MomentFitting::f_x(position2, order2-1, 0.1, 0.3)*point2[1];
                    }
                }
                BOOST_CHECK_LT(std::abs(numerical_integral), 1e-12);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(MomentFittingTestLegendrePolynomials) {
    std::cout << "Testing :: Test Moment Fitting :: Legendre Polynomials" << std::endl;
    for(int order = 1; order <= 9; ++order){

        auto ips = IntegrationPointFactory::GetIntegrationPoints(order, IntegrationPointFactory::IntegrationMethod::Gauss);
        double numerical_integral = 0.0;
        for( auto& point : ips){
            double position = point[0]* 0.2 + 0.1;
            numerical_integral += MomentFitting::f_x(position, order-1, 0.1, 0.3)*point[1]*0.2;
        }
        double analytical_int = MomentFitting::f_x_integral(0.3,order-1, 0.1, 0.3) - MomentFitting::f_x_integral(0.1,order-1, 0.1, 0.3);
        BOOST_CHECK_LT(std::abs(analytical_int-numerical_integral), 1e-12);
    }

}
BOOST_AUTO_TEST_CASE(MomentFittingSurfaceIntegral) {
    std::cout << "Testing :: Test Moment Fitting :: Surface Integral p=2" << std::endl;

    std::array<double,3> point_a_outer = {0.0, -0.1, -0.1};
    std::array<double,3> point_b_outer = {2.1, 2.1, 3.1};

    auto geometry_outer = CubeModeler::make_cube_3(point_a_outer, point_b_outer);

    std::array<double,3> point_a_inner = {-0.1, 0.0, 0.0};
    std::array<double,3> point_b_inner = {2.0, 2.0, 3.0};

    auto geometry_inner = CubeModeler::make_cube_3(point_a_inner, point_b_inner);

    std::array<double, 3> point_A = {0.0, 0.0, 0.0};
    std::array<double, 3> point_B = {2.0, 2.0, 3.0};
    std::array<int, 3> number_of_elements = {1, 1, 1};
    std::array<int, 3> order = {2, 2, 2};

    int point_distribution_factor = 3;
    double initial_triangle_edge_length = 1;
    int minimum_number_of_triangles = 20000;
    double moment_fitting_residual = 1e-8;
    std::string integration_method = "Gauss";
    int echo_level = 0;

    Parameters param(point_A, point_B, number_of_elements, order, initial_triangle_edge_length,
        minimum_number_of_triangles, moment_fitting_residual, point_distribution_factor, integration_method, echo_level);
    param.SetUseCustomizedTrimmedPointsPositionFlag(true);

    point_B = {1.0, 1.0, 1.0};
    Element element(1, point_A, point_B, param);
    auto intersection_mesh = MeshUtilities::GetIntersection(geometry_outer, geometry_inner, param);
    element.pSetSurfaceMesh(intersection_mesh);

    //auto status = EmbeddingUtilities::ComputeIntersectionMesh( *geometry_outer, *geometry_inner, element, param);

    element.GetIntegrationPointsTrimmed().clear();

    IntegrationPointUtilities::CreateGaussLegendrePoints(element.GetIntegrationPointsTrimmed(),
        point_A, point_B, order[0], order[1], order[2]);

    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPointsTrimmed() ){
        point.SetWeight(0.0);
    }

    IntegrationPointUtilities::CreateGaussLegendrePoints(element.GetIntegrationPointsInside(),
        point_A, point_B, order[0], order[1], order[2]);

    // Run Moment Fitting
    MomentFitting::ComputeReducedPointsSurfaceIntegral(element, point_distribution_factor, param);

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
    std::cout << "Error Norm: " << 1.0/27.0*std::sqrt(error_norm) << std::endl;

    BOOST_CHECK_LT(1.0/27.0*std::sqrt(error_norm), 1e-8);
} // End Testcase


BOOST_AUTO_TEST_CASE(MomentFittingSurfaceIntegralP3) {
    std::cout << "Testing :: Test Moment Fitting :: Surface Integral p=3" << std::endl;

    std::array<double,3> point_a_outer = {0.0, -0.1, -0.1};
    std::array<double,3> point_b_outer = {1.1, 1.1, 1.1};

    auto geometry_outer = CubeModeler::make_cube_3(point_a_outer, point_b_outer);

    std::array<double,3> point_a_inner = {-0.1, 0.0, 0.0};
    std::array<double,3> point_b_inner = {1.0, 1.0, 1.0};

    auto geometry_inner = CubeModeler::make_cube_3(point_a_inner, point_b_inner);

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
    auto intersection_mesh  = MeshUtilities::GetIntersection(geometry_outer, geometry_inner, param);
    element.pSetSurfaceMesh(intersection_mesh);

    element.GetIntegrationPointsTrimmed().clear();

    IntegrationPointUtilities::CreateGaussLegendrePoints(element.GetIntegrationPointsTrimmed(),
        point_A, point_B, order[0], order[1], order[2]);

    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPointsTrimmed() ){
        point.SetWeight(0.0);
    }

    IntegrationPointUtilities::CreateGaussLegendrePoints(element.GetIntegrationPointsInside(),
        point_A, point_B, order[0], order[1], order[2]);

    // Run Moment Fitting
    MomentFitting::ComputeReducedPointsSurfaceIntegral(element, point_distribution_factor, param);

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
    std::cout << "Error Norm: " << 1.0/64.0*std::sqrt(error_norm) << std::endl;

    BOOST_CHECK_LT(1.0/64.0*std::sqrt(error_norm), 1e-7);
} // End Testcase

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

//     IntegrationPointUtilities::CreateGaussLegendrePoints(element.GetIntegrationPointsTrimmedReduced(),
//         point_A, point_B, order[0], order[1], order[2]);

//     // Make sure weights are disturbed.
//     for( auto& point : element.GetIntegrationPointsTrimmedReduced() ){
//         point.SetWeight(0.0);
//     }

//     IntegrationPointUtilities::CreateGaussLegendrePoints(element.GetIntegrationPointsInside(),
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

} // End namespace Testing