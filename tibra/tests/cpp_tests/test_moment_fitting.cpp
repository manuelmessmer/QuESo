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
#include "utilities/mesh_utilities.h"
#include "io/io_utilities.h"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( MomentFittingTestSuite )

BOOST_AUTO_TEST_CASE(MomentFittingP2) {
    TIBRA_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=2" << std::endl;
    typedef boost::numeric::ublas::vector<double> VectorType;

    Parameters parameters( {Component("lower_bound", PointType(0.0, 0.0, 0.0)),
                            Component("upper_bound", PointType(2.0, 2.0, 3.0)),
                            Component("number_of_elements", Vector3i(1, 1, 1)),
                            Component("polynomial_order", Vector3i(2, 2, 2)),
                            Component("moment_fitting_residual", 1e-8),
                            Component("min_num_boundary_triangles", 5UL),
                            Component("init_point_distribution_factor", 3UL),
                            Component("integration_method", IntegrationMethod::Gauss),
                            Component("use_customized_trimmed_points", false) });

    Element element(1, {0, 0, 0}, {0.5, 0.5, 1}, parameters);

    // Construct cube over domian.
    PointType point_a_domain = {0.0, 0.0, 0.0};
    PointType point_b_domain = {1.0,1.0, 3.0};
    auto p_triangle_mesh = MeshUtilities::pGetCuboid(point_a_domain, point_b_domain);

    //MeshUtilities::Refine(*p_triangle_mesh, parameters.MinimumNumberOfTriangles() );
    auto p_boundary_ips = MakeUnique<TrimmedDomainBase::BoundaryIPVectorType>();
    for( IndexType triangle_id = 0; triangle_id < p_triangle_mesh->NumOfTriangles(); ++triangle_id ) {
            IndexType method = 3; // This will create 6 points per triangle.
            auto p_new_points = p_triangle_mesh->pGetIPsGlobal(triangle_id, method);
            p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
    }

    // Distribtue Gauss points within element.
    element.GetIntegrationPoints().clear();
    SingleElement::AssembleIPs(element, parameters);
    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    IO::WriteMeshToSTL(*p_triangle_mesh, "mesh.stl", true);
    // Run Moment Fitting
    VectorType constant_terms{};
    MomentFitting::ComputeConstantTerms(element, p_boundary_ips, constant_terms, parameters);
    MomentFitting::MomentFitting1(constant_terms, element.GetIntegrationPoints(), element, parameters);
    auto& points_moment_fitting = element.GetIntegrationPoints();

    // Get Gauss points as reference
    Element::IntegrationPointVectorType points_gauss_legendre{};
    SingleElement::AssembleIPs(points_gauss_legendre, element.GetLowerBoundParam(),
        element.GetUpperBoundParam(), {2, 2, 2});

    double error_norm = 0.0;
    // Check if weights are similar
    for( int i = 0; i < 27; ++i){
        double weight_mf = points_moment_fitting[i].GetWeight();
        double weight_gl = points_gauss_legendre[i].GetWeight();
        double error = (weight_mf-weight_gl)/weight_gl;
        error_norm += std::pow(error,2);
        BOOST_CHECK_CLOSE_FRACTION(weight_mf, weight_gl, 1e-12);
    }

    BOOST_CHECK_LT(1.0/27.0*std::sqrt(error_norm), 1e-10);
} // End Testcase


BOOST_AUTO_TEST_CASE(MomentFittingP3) {
    TIBRA_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=3" << std::endl;
    typedef boost::numeric::ublas::vector<double> VectorType;

    Parameters parameters( {Component("lower_bound", PointType(0.0, 0.0, 0.0)),
                            Component("upper_bound", PointType(2.0, 2.0, 1.0)),
                            Component("number_of_elements", Vector3i(1, 1, 1)),
                            Component("polynomial_order", Vector3i(3, 3, 3)),
                            Component("moment_fitting_residual", 1e-8),
                            Component("min_num_boundary_triangles", 500UL),
                            Component("init_point_distribution_factor", 3UL),
                            Component("integration_method", IntegrationMethod::Gauss),
                            Component("use_customized_trimmed_points", false) });

    Element element(1, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, parameters);

    // Construct cube over domian.
    PointType point_a_domain = {0.0, 0.0, 0.0};
    PointType point_b_domain = {2.0, 2.0, 1.0};
    auto p_triangle_mesh = MeshUtilities::pGetCuboid(point_a_domain, point_b_domain);

    MeshUtilities::Refine(*p_triangle_mesh, parameters.MinimumNumberOfTriangles() );
    auto p_boundary_ips = MakeUnique<TrimmedDomainBase::BoundaryIPVectorType>();
    for( IndexType triangle_id = 0; triangle_id < p_triangle_mesh->NumOfTriangles(); ++triangle_id ) {
            IndexType method = 3; // This will create 6 points per triangle.
            auto p_new_points = p_triangle_mesh->pGetIPsGlobal(triangle_id, method);
            p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
    }

    // Distribtue Gauss points within element.
    element.GetIntegrationPoints().clear();
    SingleElement::AssembleIPs(element, parameters);
    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    IO::WriteMeshToSTL(*p_triangle_mesh, "mesh.stl", true);
    // Run Moment Fitting
    VectorType constant_terms{};
    MomentFitting::ComputeConstantTerms(element, p_boundary_ips, constant_terms, parameters);
    MomentFitting::MomentFitting1(constant_terms, element.GetIntegrationPoints(), element, parameters);
    auto& points_moment_fitting = element.GetIntegrationPoints();

    // Get Gauss points as reference
    Element::IntegrationPointVectorType points_gauss_legendre{};
    SingleElement::AssembleIPs(points_gauss_legendre, element.GetLowerBoundParam(),
        element.GetUpperBoundParam(), {3, 3, 3});

    double error_norm = 0.0;
    // Check if weights are similar
    for( int i = 0; i < 64; ++i){
        double weight_mf = points_moment_fitting[i].GetWeight();
        double weight_gl = points_gauss_legendre[i].GetWeight();
        double error = (weight_mf-weight_gl)/weight_gl;
        error_norm += std::pow(error,2);
        BOOST_CHECK_CLOSE_FRACTION(weight_mf, weight_gl, 1e-6);
    }
    //TIBRA_INFO << "Error Norm: " << 1.0/64.0*std::sqrt(error_norm) << std::endl;

    BOOST_CHECK_LT(1.0/64.0*std::sqrt(error_norm), 1e-7);
} // End Testcase

BOOST_AUTO_TEST_CASE(MomentFittingP4) {
    TIBRA_INFO << "Testing :: Test Moment Fitting :: Surface Integral p=3" << std::endl;
    typedef boost::numeric::ublas::vector<double> VectorType;

    Parameters parameters( {Component("lower_bound", PointType(0.0, 0.0, 0.0)),
                            Component("upper_bound", PointType(2.0, 2.0, 1.0)),
                            Component("number_of_elements", Vector3i(1, 1, 1)),
                            Component("polynomial_order", Vector3i(4, 4, 4)),
                            Component("moment_fitting_residual", 1e-8),
                            Component("min_num_boundary_triangles", 1000UL),
                            Component("init_point_distribution_factor", 3UL),
                            Component("integration_method", IntegrationMethod::Gauss),
                            Component("use_customized_trimmed_points", false) });

    Element element(1, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, parameters);

    // Construct cube over domian.
    PointType point_a_domain = {0.0, 0.0, 0.0};
    PointType point_b_domain = {2.0, 2.0, 1.0};
    auto p_triangle_mesh = MeshUtilities::pGetCuboid(point_a_domain, point_b_domain);

    MeshUtilities::Refine(*p_triangle_mesh, parameters.MinimumNumberOfTriangles() );
    auto p_boundary_ips = MakeUnique<TrimmedDomainBase::BoundaryIPVectorType>();
    for( IndexType triangle_id = 0; triangle_id < p_triangle_mesh->NumOfTriangles(); ++triangle_id ) {
            IndexType method = 3; // This will create 6 points per triangle.
            auto p_new_points = p_triangle_mesh->pGetIPsGlobal(triangle_id, method);
            p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
    }

    // Distribtue Gauss points within element.
    element.GetIntegrationPoints().clear();
    SingleElement::AssembleIPs(element, parameters);
    // Make sure weights are disturbed.
    for( auto& point : element.GetIntegrationPoints() ){
        point.SetWeight(0.0);
    }

    IO::WriteMeshToSTL(*p_triangle_mesh, "mesh.stl", true);
    // Run Moment Fitting
    VectorType constant_terms{};
    MomentFitting::ComputeConstantTerms(element, p_boundary_ips, constant_terms, parameters);
    MomentFitting::MomentFitting1(constant_terms, element.GetIntegrationPoints(), element, parameters);
    auto& points_moment_fitting = element.GetIntegrationPoints();

    // Get Gauss points as reference
    Element::IntegrationPointVectorType points_gauss_legendre{};
    SingleElement::AssembleIPs(points_gauss_legendre, element.GetLowerBoundParam(),
        element.GetUpperBoundParam(), {4, 4, 4});

    double error_norm = 0.0;
    // Check if weights are similar
    for( int i = 0; i < 125; ++i){
        double weight_mf = points_moment_fitting[i].GetWeight();
        double weight_gl = points_gauss_legendre[i].GetWeight();
        double error = (weight_mf-weight_gl)/weight_gl;
        error_norm += std::pow(error,2);
        BOOST_CHECK_CLOSE_FRACTION(weight_mf, weight_gl, 1e-6);
    }
    //TIBRA_INFO << "Error Norm: " << 1.0/64.0*std::sqrt(error_norm) << std::endl;

    BOOST_CHECK_LT(1.0/125.0*std::sqrt(error_norm), 1e-7);
} // End Testcase

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra
