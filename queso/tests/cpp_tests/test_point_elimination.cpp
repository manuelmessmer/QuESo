// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "includes/checks.hpp"
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

    Parameters parameters( {Component("lower_bound_xyz", Vector3d(-1.5, -1.5, -1.0)),
                            Component("upper_bound_xyz", Vector3d(1.5, 1.5, 12.0)),
                            Component("lower_bound_uvw", Vector3d(0.0, 0.0, 0.0)),
                            Component("upper_bound_uvw", Vector3d(1.0, 1.0, 1.0)),
                            Component("number_of_elements", Vector3i(6, 6, 13))} );

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cylinder.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 1e-3;
    const IndexType min_num_triangles = 500;

    Mapper mapper(parameters);
    IndexType number_trimmed_elements = 0;
    for( IndexType i = 0; i < mapper.NumberOfElements(); ++i){
        const BoundingBoxType bounding_box = mapper.GetBoundingBoxXYZFromIndex(i);
        Vector3d lower_bound_xyz = bounding_box.first;
        Vector3d upper_bound_xyz = bounding_box.second;

        const BoundingBoxType bounding_box_uvw = mapper.GetBoundingBoxUVWFromIndex(i);
        Vector3d lower_bound_uvw = bounding_box_uvw.first;
        Vector3d upper_bound_uvw = bounding_box_uvw.second;

        // Construct element
        Element element(1, MakeBox(lower_bound_xyz, upper_bound_xyz),
                           MakeBox(lower_bound_uvw, upper_bound_uvw));

        if( brep_operator.GetIntersectionState(lower_bound_xyz, upper_bound_xyz) == IntersectionStatus::Trimmed){
            // Get trimmed domain
            auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound_xyz, upper_bound_xyz, min_vol_ratio, min_num_triangles);
            if( p_trimmed_domain ){
                ++number_trimmed_elements;

                // Add trimmed domain to element.
                element.SetIsTrimmed(true);
                element.pSetTrimmedDomain(p_trimmed_domain);

                // Run point elimination
                const auto residual = QuadratureTrimmedElementTester::AssembleIPs(element, rOrder, Residual);

                // Check if residual is smaller than targeted.
                QuESo_CHECK_LT(residual, 1e-6);

                // Must be more points than p*p*p.
                auto& r_points = element.GetIntegrationPoints();

                QuESo_CHECK_LT(r_points.size(), (rOrder[0]+1)*(rOrder[1]+1)*(rOrder[2]+1)+1);
                QuESo_CHECK_GT(r_points.size(), rOrder[0]*rOrder[1]*rOrder[2]);

                // Get copy of points.
                Element::IntegrationPointVectorType copy_points(r_points);

                // Compute constant terms.
                std::vector<double> constant_terms{};
                auto p_boundary_ips = element.pGetTrimmedDomain()->pGetBoundaryIps();
                QuadratureTrimmedElementTester::ComputeConstantTerms(constant_terms, p_boundary_ips, element, rOrder);

                // Run moment fitting again.
                const auto residual_2 = QuadratureTrimmedElementTester::MomentFitting(constant_terms, r_points, element, rOrder);

                // Check if residual and weights are the same.
                QuESo_CHECK_NEAR( residual, residual_2, EPS4 );
                double volume = 0.0;
                for( IndexType i = 0; i < r_points.size(); ++i){
                    const double weight1 = r_points[i].Weight();
                    QuESo_CHECK_GT(weight1, EPS4);
                    const double weight2 = copy_points[i].Weight();
                    QuESo_CHECK_RELATIVE_NEAR( weight1, weight2 , EPS2 );
                    volume += weight1*element.DetJ(); // Multiplied with det(J).
                }
                // Check if integration points contain correct volume;
                const auto& r_mesh = element.pGetTrimmedDomain()->GetTriangleMesh();
                const double ref_volume = MeshUtilities::Volume(r_mesh);
                QuESo_CHECK_RELATIVE_NEAR(volume, ref_volume, 100.0*Residual);
            }
        }
    }
    QuESo_CHECK_EQUAL(number_trimmed_elements, 120);
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

    Parameters parameters( {Component("lower_bound_xyz", Vector3d(-130.0, -110.0, -110.0)),
                            Component("upper_bound_xyz", Vector3d(-40, 10.0, 10.0)),
                            Component("lower_bound_uvw", Vector3d(-130.0, -110.0, -110.0)),
                            Component("upper_bound_uvw", Vector3d(-40, 10.0, 10.0)),
                            Component("number_of_elements", Vector3i(9, 12, 12)) } );


    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/steering_knuckle.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 1e-3;
    const IndexType min_num_triangles = 500;

    Mapper mapper(parameters);
    IndexType number_trimmed_elements = 0;
    for( IndexType i = 0; i < mapper.NumberOfElements(); ++i){
        const BoundingBoxType bounding_box = mapper.GetBoundingBoxXYZFromIndex(i);
        Vector3d lower_bound_xyz = bounding_box.first;
        Vector3d upper_bound_xyz = bounding_box.second;

        const BoundingBoxType bounding_box_uvw = mapper.GetBoundingBoxUVWFromIndex(i);
        Vector3d lower_bound_uvw = bounding_box_uvw.first;
        Vector3d upper_bound_uvw = bounding_box_uvw.second;

        // Construct element
        Element element(1, MakeBox(lower_bound_xyz, upper_bound_xyz),
                           MakeBox(lower_bound_uvw, upper_bound_uvw));

        if( brep_operator.GetIntersectionState(lower_bound_xyz, upper_bound_xyz) == IntersectionStatus::Trimmed){
            // Get trimmed domain
            auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound_xyz, upper_bound_xyz, min_vol_ratio, min_num_triangles);
            if( p_trimmed_domain ){
                ++number_trimmed_elements;

                // Add trimmed domain to element.
                element.SetIsTrimmed(true);
                element.pSetTrimmedDomain(p_trimmed_domain);

                // Run point elimination
                const auto residual = QuadratureTrimmedElementTester::AssembleIPs(element, {2, 2, 2}, 1e-8);

                // Check if residual is smaller than targeted.
                QuESo_CHECK_LT(residual, 1e-8);

                // Must be more points than p*p*p.
                auto& r_points = element.GetIntegrationPoints();
                QuESo_CHECK_LT(r_points.size(), 28);
                QuESo_CHECK_GT(r_points.size(), 7);

                // Get copy of points.
                Element::IntegrationPointVectorType copy_points(r_points);

                // Compute constant terms.
                std::vector<double> constant_terms{};
                auto p_boundary_ips = element.pGetTrimmedDomain()->pGetBoundaryIps();
                QuadratureTrimmedElementTester::ComputeConstantTerms(constant_terms, p_boundary_ips, element, {2, 2, 2});

                // Run moment fitting again.
                const auto residual_2 = QuadratureTrimmedElementTester::MomentFitting(constant_terms, r_points, element, {2, 2, 2});

                // Check if residual and weights are the same.
                QuESo_CHECK_LT( residual, residual_2+EPS4 );
                QuESo_CHECK_LT( residual_2, residual+EPS4 );
                double volume = 0.0;
                for( IndexType i = 0; i < r_points.size(); ++i){
                    const double weight1 = r_points[i].Weight();
                    QuESo_CHECK_GT(weight1, EPS4);
                    const double weight2 = copy_points[i].Weight();
                    const double error = std::abs(weight1 - weight2)/ weight1;
                    QuESo_CHECK_LT( error , EPS2 );
                    volume += weight1*element.DetJ(); // Multiplied with det(J).
                }
                // Check if integration points contain correct volume;
                const auto& r_mesh = element.pGetTrimmedDomain()->GetTriangleMesh();
                const double ref_volume = MeshUtilities::Volume(r_mesh);
                const double volume_error = std::abs(volume - ref_volume)/ ref_volume;
                QuESo_CHECK_LT(volume_error, 1e-6); // Note can not be better as moment fitting residual.
            }
        }
    }
    QuESo_CHECK_EQUAL(number_trimmed_elements, 80);
}

BOOST_AUTO_TEST_CASE(PointEliminationElephantTest) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Elephant" << std::endl;

    Parameters parameters( {Component("lower_bound_xyz", Vector3d(-0.4, -0.6, -0.35)),
                            Component("upper_bound_xyz", Vector3d(0.4, 0.6, 0.35)),
                            Component("lower_bound_uvw", Vector3d(-1.0, -1.0, -1.0)),
                            Component("upper_bound_uvw", Vector3d( 1.0,  1.0,  1.0)),
                            Component("b_spline_mesh", false),
                            Component("number_of_elements", Vector3i(8, 12, 7)) } );

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/elephant.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 1e-3;
    const IndexType min_num_triangles = 500;

    Mapper mapper(parameters);
    IndexType number_trimmed_elements = 0;
    for( IndexType i = 0; i < mapper.NumberOfElements(); ++i){
        const BoundingBoxType bounding_box = mapper.GetBoundingBoxXYZFromIndex(i);
        Vector3d lower_bound_xyz = bounding_box.first;
        Vector3d upper_bound_xyz = bounding_box.second;

        const BoundingBoxType bounding_box_uvw = mapper.GetBoundingBoxUVWFromIndex(i);
        Vector3d lower_bound_uvw = bounding_box_uvw.first;
        Vector3d upper_bound_uvw = bounding_box_uvw.second;

        // Construct element
        Element element(1, MakeBox(lower_bound_xyz, upper_bound_xyz),
                           MakeBox(lower_bound_uvw, upper_bound_uvw));

        if( brep_operator.GetIntersectionState(lower_bound_xyz, upper_bound_xyz) == IntersectionStatus::Trimmed){
            // Get trimmed domain
            auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound_xyz, upper_bound_xyz, min_vol_ratio, min_num_triangles);
            if( p_trimmed_domain ){
                ++number_trimmed_elements;

                // Add trimmed domain to element.
                element.SetIsTrimmed(true);
                element.pSetTrimmedDomain(p_trimmed_domain);

                // Run point elimination
                const auto residual = QuadratureTrimmedElementTester::AssembleIPs(element, {2, 2, 2}, 1e-8);

                // Check if residual is smaller than targeted.
                QuESo_CHECK_LT(residual, 1e-8);

                // Must be more points than p*p*p.
                auto& r_points = element.GetIntegrationPoints();
                QuESo_CHECK_LT(r_points.size(), 28);
                QuESo_CHECK_GT(r_points.size(), 7);

                // Get copy of points.
                Element::IntegrationPointVectorType copy_points(r_points);

                // Compute constant terms.
                std::vector<double> constant_terms{};
                auto p_boundary_ips = element.pGetTrimmedDomain()->pGetBoundaryIps();
                QuadratureTrimmedElementTester::ComputeConstantTerms(constant_terms, p_boundary_ips, element, {2, 2, 2});

                // Run moment fitting again.
                const auto residual_2 = QuadratureTrimmedElementTester::MomentFitting(constant_terms, r_points, element, {2, 2, 2});

                // Check if residual and weights are the same.
                QuESo_CHECK_LT( residual, residual_2+EPS4 );
                QuESo_CHECK_LT( residual_2, residual+EPS4 );
                double volume = 0.0;
                for( IndexType i = 0; i < r_points.size(); ++i){
                    const double weight1 = r_points[i].Weight();
                    QuESo_CHECK_GT(weight1, EPS4);
                    const double weight2 = copy_points[i].Weight();
                    const double error = std::abs(weight1 - weight2)/ weight1;
                    QuESo_CHECK_LT( error , EPS2 );
                    volume += weight1*element.DetJ(); // Multiplied with det(J).
                }

                // Check if integration points contain correct volume;
                const auto& r_mesh = element.pGetTrimmedDomain()->GetTriangleMesh();
                const double ref_volume = MeshUtilities::Volume(r_mesh);
                const double volume_error = std::abs(volume - ref_volume);
                QuESo_CHECK_LT(volume_error, 1e-7);
            }
        }
    }
    QuESo_CHECK_EQUAL(number_trimmed_elements, 153);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso