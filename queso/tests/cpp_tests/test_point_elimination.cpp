// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "containers/element_container.hpp"
#include "quadrature/trimmed_element.h"
#include "containers/triangle_mesh.hpp"
#include "embedding/brep_operator.h"
#include "io/io_utilities.h"
#include "tests/cpp_tests/class_testers/trimmed_element_tester.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( PointEliminationTestSuite )

void RunCylinder(const Vector3i& rOrder, double Residual){
    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};
    Vector3d lower_bound = {-1.5, -1.5, -1.0 };
    Vector3d upper_bound = {1.5, 1.5, 12.0 };
    Parameters parameters( {Component("lower_bound_xyz", lower_bound),
                            Component("upper_bound_xyz", upper_bound),
                            Component("lower_bound_uvw", lower_bound),
                            Component("upper_bound_uvw", upper_bound),
                            Component("min_num_boundary_triangles", 500UL),
                            Component("moment_fitting_residual", Residual),
                            Component("number_of_elements", number_of_elements),
                            Component("init_point_distribution_factor", 1UL),
                            Component("echo_level", 1UL),
                            Component("min_element_volume_ratio", 1e-3),
                            Component("polynomial_order", rOrder ) } );

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cylinder.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh,parameters);

    const double delta_x = 0.5;
    const double delta_y = 0.5;
    const double delta_z = 1;

    IndexType number_trimmed_elements = 0;
    for(double x = lower_bound[0]; x <= upper_bound[0]; x += delta_x){
        for(double y = lower_bound[1]; y <= upper_bound[1]; y += delta_y){
            for(double z = lower_bound[2]; z <= upper_bound[2]; z += delta_z){

                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                // Construct element
                Element element(1, MakeBox(local_lower_bound, local_upper_bound),
                                   MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), parameters);

                if( brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound) == IntersectionStatus::Trimmed){
                    // Get trimmed domain
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(local_lower_bound, local_upper_bound);
                    if( p_trimmed_domain ){
                        ++number_trimmed_elements;

                        // Add trimmed domain to element.
                        element.SetIsTrimmed(true);
                        element.pSetTrimmedDomain(p_trimmed_domain);

                        // Run point elimination
                        const auto residual = QuadratureTrimmedElementTester::AssembleIPs(element, parameters);

                        // Check if residual is smaller than targeted.
                        BOOST_CHECK_LT(residual, 1e-6);

                        // Must be more points than p*p*p.
                        auto& r_points = element.GetIntegrationPoints();

                        BOOST_CHECK_LT(r_points.size(), (rOrder[0]+1)*(rOrder[1]+1)*(rOrder[2]+1)+1);
                        BOOST_CHECK_GT(r_points.size(), rOrder[0]*rOrder[1]*rOrder[2]);

                        // Get copy of points.
                        Element::IntegrationPointVectorType copy_points(r_points);

                        // Compute constant terms.
                        VectorType constant_terms{};
                        auto p_boundary_ips = element.pGetTrimmedDomain()->pGetBoundaryIps();
                        QuadratureTrimmedElementTester::ComputeConstantTerms(constant_terms, p_boundary_ips, element, parameters);

                        // Run moment fitting again.
                        const auto residual_2 = QuadratureTrimmedElementTester::MomentFitting(constant_terms, r_points, element, parameters);

                        // Check if residual and weights are the same.
                        BOOST_CHECK_LT( residual, residual_2+EPS4 );
                        BOOST_CHECK_LT( residual_2, residual+EPS4 );
                        double volume = 0.0;
                        for( IndexType i = 0; i < r_points.size(); ++i){
                            const double weight1 = r_points[i].GetWeight();
                            BOOST_CHECK_GT(weight1, EPS4);
                            const double weight2 = copy_points[i].GetWeight();
                            const double error = std::abs(weight1 - weight2)/ weight1;
                            BOOST_CHECK_LT( error , EPS2 );
                            volume += weight1*element.DetJ(); // Multiplied with det(J).
                        }
                        // Check if integration points contain correct volume;
                        const auto& r_mesh = element.pGetTrimmedDomain()->GetTriangleMesh();
                        const double ref_volume = MeshUtilities::Volume(r_mesh);
                        const double volume_error = std::abs(volume - ref_volume)/ ref_volume;
                        BOOST_CHECK_LT(volume_error, Residual*100.0); // Note can not be better as moment fitting residual.
                    }
                }
            }
        }
    }
    BOOST_CHECK_EQUAL(number_trimmed_elements, 120);
}

BOOST_AUTO_TEST_CASE(PointEliminationCylinder1Test) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Cylinder Quadratic" << std::endl;
    RunCylinder({2, 2, 2}, 1e-8);
}

BOOST_AUTO_TEST_CASE(PointEliminationCylinder2Test) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Cylinder Cubic" << std::endl;
    RunCylinder({3, 3, 3}, 1e-8);
}

BOOST_AUTO_TEST_CASE(PointEliminationCylinder4Test) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Cylinder Mixed" << std::endl;
    RunCylinder({2, 3, 4}, 1e-7);
}

BOOST_AUTO_TEST_CASE(PointEliminationKnuckleTest) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Knuckle" << std::endl;

    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};
    Vector3d lower_bound = {-130.0, -110.0, -110.0 };
    Vector3d upper_bound = {-50, 0.0, 0.0 };
    Parameters parameters( {Component("lower_bound_xyz", lower_bound),
                            Component("upper_bound_xyz", upper_bound),
                            Component("lower_bound_uvw", lower_bound),
                            Component("upper_bound_uvw", upper_bound),
                            Component("min_num_boundary_triangles", 500UL),
                            Component("moment_fitting_residual", 1e-8),
                            Component("number_of_elements", number_of_elements),
                            Component("init_point_distribution_factor", 1UL),
                            Component("echo_level", 1UL),
                            Component("min_element_volume_ratio", 1e-3),
                            Component("polynomial_order", Vector3i(2,2,2) ) } );


    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/steering_knuckle.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh,parameters);

    const double delta_x = 10;
    const double delta_y = 10;
    const double delta_z = 10;

    IndexType number_trimmed_elements = 0;
    for(double x = lower_bound[0]; x <= upper_bound[0]; x += delta_x){
        for(double y = lower_bound[1]; y <= upper_bound[1]; y += delta_y){
            for(double z = lower_bound[2]; z <= upper_bound[2]; z += delta_z){

                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                // Construct element
                Element element(1, MakeBox(local_lower_bound, local_upper_bound),
                                   MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), parameters);

                if( brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound) == IntersectionStatus::Trimmed){
                    // Get trimmed domain
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(local_lower_bound, local_upper_bound);
                    if( p_trimmed_domain ){
                        ++number_trimmed_elements;

                        // Add trimmed domain to element.
                        element.SetIsTrimmed(true);
                        element.pSetTrimmedDomain(p_trimmed_domain);

                        // Run point elimination
                        const auto residual = QuadratureTrimmedElementTester::AssembleIPs(element, parameters);

                        // Check if residual is smaller than targeted.
                        BOOST_CHECK_LT(residual, 1e-8);

                        // Must be more points than p*p*p.
                        auto& r_points = element.GetIntegrationPoints();
                        BOOST_CHECK_LT(r_points.size(), 28);
                        BOOST_CHECK_GT(r_points.size(), 7);

                        // Get copy of points.
                        Element::IntegrationPointVectorType copy_points(r_points);

                        // Compute constant terms.
                        VectorType constant_terms{};
                        auto p_boundary_ips = element.pGetTrimmedDomain()->pGetBoundaryIps();
                        QuadratureTrimmedElementTester::ComputeConstantTerms(constant_terms, p_boundary_ips, element, parameters);

                        // Run moment fitting again.
                        const auto residual_2 = QuadratureTrimmedElementTester::MomentFitting(constant_terms, r_points, element, parameters);

                        // Check if residual and weights are the same.
                        BOOST_CHECK_LT( residual, residual_2+EPS4 );
                        BOOST_CHECK_LT( residual_2, residual+EPS4 );
                        double volume = 0.0;
                        for( IndexType i = 0; i < r_points.size(); ++i){
                            const double weight1 = r_points[i].GetWeight();
                            BOOST_CHECK_GT(weight1, EPS4);
                            const double weight2 = copy_points[i].GetWeight();
                            const double error = std::abs(weight1 - weight2)/ weight1;
                            BOOST_CHECK_LT( error , EPS2 );
                            volume += weight1*element.DetJ(); // Multiplied with det(J).
                        }
                        // Check if integration points contain correct volume;
                        const auto& r_mesh = element.pGetTrimmedDomain()->GetTriangleMesh();
                        const double ref_volume = MeshUtilities::Volume(r_mesh);
                        const double volume_error = std::abs(volume - ref_volume)/ ref_volume;
                        BOOST_CHECK_LT(volume_error, 1e-6); // Note can not be better as moment fitting residual.
                    }
                }
            }
        }
    }
    BOOST_CHECK_EQUAL(number_trimmed_elements, 80);
}

BOOST_AUTO_TEST_CASE(PointEliminationElephantTest) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Elephant" << std::endl;

    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};
    Vector3d lower_bound = {-0.4, -0.6, -0.35 };
    Vector3d upper_bound = {0.4, 0.6, 0.35 };
    Parameters parameters( {Component("lower_bound_xyz", lower_bound),
                            Component("upper_bound_xyz", upper_bound),
                            Component("lower_bound_uvw", lower_bound),
                            Component("upper_bound_uvw", upper_bound),
                            Component("min_num_boundary_triangles", 500UL),
                            Component("moment_fitting_residual", 1e-8),
                            Component("number_of_elements", number_of_elements),
                            Component("init_point_distribution_factor", 1UL),
                            Component("echo_level", 1UL),
                            Component("min_element_volume_ratio", 1e-3),
                            Component("polynomial_order", Vector3i(2,2,2) ) } );

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/elephant.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh,parameters);

    const double delta_x = 0.1;
    const double delta_y = 0.1;
    const double delta_z = 0.1;

    IndexType number_trimmed_elements = 0;
    for(double x = lower_bound[0]; x <= upper_bound[0]; x += delta_x){
        for(double y = lower_bound[1]; y <= upper_bound[1]; y += delta_y){
            for(double z = lower_bound[2]; z <= upper_bound[2]; z += delta_z){

                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                // Construct element
                Element element(1, MakeBox(local_lower_bound, local_upper_bound),
                                   MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), parameters);

                if( brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound) == IntersectionStatus::Trimmed){
                    // Get trimmed domain
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(local_lower_bound, local_upper_bound);
                    if( p_trimmed_domain ){
                        ++number_trimmed_elements;

                        // Add trimmed domain to element.
                        element.SetIsTrimmed(true);
                        element.pSetTrimmedDomain(p_trimmed_domain);

                        // Run point elimination
                        const auto residual = QuadratureTrimmedElementTester::AssembleIPs(element, parameters);

                        // Check if residual is smaller than targeted.
                        BOOST_CHECK_LT(residual, 1e-8);

                        // Must be more points than p*p*p.
                        auto& r_points = element.GetIntegrationPoints();
                        BOOST_CHECK_LT(r_points.size(), 28);
                        BOOST_CHECK_GT(r_points.size(), 7);

                        // Get copy of points.
                        Element::IntegrationPointVectorType copy_points(r_points);

                        // Compute constant terms.
                        VectorType constant_terms{};
                        auto p_boundary_ips = element.pGetTrimmedDomain()->pGetBoundaryIps();
                        QuadratureTrimmedElementTester::ComputeConstantTerms(constant_terms, p_boundary_ips, element, parameters);

                        // Run moment fitting again.
                        const auto residual_2 = QuadratureTrimmedElementTester::MomentFitting(constant_terms, r_points, element, parameters);

                        // Check if residual and weights are the same.
                        BOOST_CHECK_LT( residual, residual_2+EPS4 );
                        BOOST_CHECK_LT( residual_2, residual+EPS4 );
                        double volume = 0.0;
                        for( IndexType i = 0; i < r_points.size(); ++i){
                            const double weight1 = r_points[i].GetWeight();
                            BOOST_CHECK_GT(weight1, EPS4);
                            const double weight2 = copy_points[i].GetWeight();
                            const double error = std::abs(weight1 - weight2)/ weight1;
                            BOOST_CHECK_LT( error , EPS2 );
                            volume += weight1*element.DetJ(); // Multiplied with det(J).
                        }
                        // Check if integration points contain correct volume;
                        const auto& r_mesh = element.pGetTrimmedDomain()->GetTriangleMesh();
                        const double ref_volume = MeshUtilities::Volume(r_mesh);
                        const double volume_error = std::abs(volume - ref_volume);
                        BOOST_CHECK_LT(volume_error, 1e-7);
                    }
                }
            }
        }
    }
    BOOST_CHECK_EQUAL(number_trimmed_elements, 153);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso