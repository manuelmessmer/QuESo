// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "includes/checks.hpp"
#include "containers/element_container.hpp"
#include "containers/triangle_mesh.hpp"
#include "embedding/brep_operator.h"
#include "io/io_utilities.h"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( PointClassifierOnTrimmedDomainTestSuite )

BOOST_AUTO_TEST_CASE(CylinderPointClassifierOnTrimmedDomainTest) {

    QuESo_INFO << "Testing :: Test Point Classifier On Trimmed Domain:: Cylinder Point Classifier" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cylinder.stl");

    Parameters params( {Component("lower_bound_xyz", PointType(0.0, 0.0, 0.0)),
                        Component("upper_bound_xyz", PointType(1.0, 1.0, 1.0)),
                        Component("lower_bound_uvw", PointType(0.0, 0.0, 0.0)),
                        Component("upper_bound_uvw", PointType(1.0, 1.0, 1.0)),
                        Component("number_of_elements", Vector3i(1, 1, 1)),
                        Component("min_element_volume_ratio", 0.0) });

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh, params);

    const double delta_x = 0.50;
    const double delta_y = 0.50;
    const double delta_z = 0.50;

    IndexType num_of_trimmed_elements = 0;
    for(double x = -1.5; x <= 1.5; x += delta_x){
        for(double y = -1.5; y <= 1.5; y += delta_y){
            for(double z = -1.0; z <= 12; z += delta_z){
                Vector3d lower_bound = {x, y, z};
                Vector3d upper_bound = {x+delta_x, y+delta_y, z+delta_z};
                if( brep_operator.GetIntersectionState(lower_bound, upper_bound) == IntersectionStatus::Trimmed){
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound);
                    const double delta_p_x = 0.1;
                    const double delta_p_y = 0.1;
                    const double delta_p_z = 0.1;
                    for(double p_x = lower_bound[0]+delta_p_x/2.0; p_x <= upper_bound[0]; p_x += delta_p_x){
                        for(double p_y = lower_bound[1]+delta_p_y/2.0; p_y <= upper_bound[1]; p_y += delta_p_y){
                            for(double p_z = lower_bound[2]+delta_p_z/2.0; p_z <= upper_bound[2]; p_z += delta_p_z){
                                // Make sure both functions produce the same result.
                                PointType test_point = {p_x, p_y, p_z};
                                bool res1 = p_trimmed_domain->IsInsideTrimmedDomain(test_point);
                                bool res2 = brep_operator.IsInside(test_point);
                                QuESo_CHECK_EQUAL(res1, res2);
                            }
                        }
                    }
                    num_of_trimmed_elements++;
                }
            }
        }
    }
    QuESo_CHECK_EQUAL(num_of_trimmed_elements, 240);
}

BOOST_AUTO_TEST_CASE(CubePointClassifierOnTrimmedDomainTest) {
    QuESo_INFO << "Testing :: Test Point Classifier On Trimmed Domain:: Cube with Cavity Point Classifier" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cube_with_cavity.stl");

    Parameters params( {Component("lower_bound_xyz", PointType(0.0, 0.0, 0.0)),
                        Component("upper_bound_xyz", PointType(1.0, 1.0, 1.0)),
                        Component("lower_bound_uvw", PointType(0.0, 0.0, 0.0)),
                        Component("upper_bound_uvw", PointType(1.0, 1.0, 1.0)),
                        Component("number_of_elements", Vector3i(1, 1, 1)),
                        Component("min_element_volume_ratio", 0.0) });

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh, params);

    const double delta_x = 0.15;
    const double delta_y = 0.15;
    const double delta_z = 0.15;

    IndexType num_of_trimmed_elements = 0;
    for(double x = -1.5001; x <= 1.5; x += delta_x){
        for(double y = -1.5001; y <= 1.5; y += delta_y){
            for(double z = -1.5001; z <= 1.5; z += delta_z){
                Vector3d lower_bound = {x, y, z};
                Vector3d upper_bound = {x+delta_x, y+delta_y, z+delta_z};
                if( brep_operator.GetIntersectionState(lower_bound, upper_bound) == IntersectionStatus::Trimmed){
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound);
                    const double delta_p_x = 0.2;
                    const double delta_p_y = 0.2;
                    const double delta_p_z = 0.2;
                    for(double p_x = lower_bound[0]+delta_p_x/2.0; p_x <= upper_bound[0]; p_x += delta_p_x){
                        for(double p_y = lower_bound[1]+delta_p_y/2.0; p_y <= upper_bound[1]; p_y += delta_p_y){
                            for(double p_z = lower_bound[2]+delta_p_z/2.0; p_z <= upper_bound[2]; p_z += delta_p_z){
                                // Make sure both functions produce the same result.
                                PointType test_point = {p_x, p_y, p_z};
                                bool res1 = p_trimmed_domain->IsInsideTrimmedDomain(test_point);
                                bool res2 = brep_operator.IsInside(test_point);
                                QuESo_CHECK_EQUAL(res1, res2);
                            }
                        }
                    }
                    num_of_trimmed_elements++;
                }
            }
        }
    }
    QuESo_CHECK_EQUAL(num_of_trimmed_elements, 3226);
}

BOOST_AUTO_TEST_CASE(ElephantPointClassifierOnTrimmedDomainTest) {
    QuESo_INFO << "Testing :: Test Point Classifier On Trimmed Domain:: Elephant Point Classifier" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/elephant.stl");

    Parameters params( {Component("lower_bound_xyz", PointType(0.0, 0.0, 0.0)),
                        Component("upper_bound_xyz", PointType(1.0, 1.0, 1.0)),
                        Component("lower_bound_uvw", PointType(0.0, 0.0, 0.0)),
                        Component("upper_bound_uvw", PointType(1.0, 1.0, 1.0)),
                        Component("number_of_elements", Vector3i(1, 1, 1)),
                        Component("min_element_volume_ratio", 0.0) });

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh, params);

    const double delta_x = 0.05;
    const double delta_y = 0.05;
    const double delta_z = 0.05;

    IndexType num_of_trimmed_elements = 0;
    for(double x = -0.4; x <= 0.4; x += delta_x){
        for(double y = -0.6; y <= 0.6; y += delta_y){
            for(double z = -0.35; z <= 0.35; z += delta_x){
                Vector3d lower_bound = {x, y, z};
                Vector3d upper_bound = {x+delta_x, y+delta_y, z+delta_z};
                if( brep_operator.GetIntersectionState(lower_bound, upper_bound) == IntersectionStatus::Trimmed){
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound);
                    const double delta_p_x = 0.01;
                    const double delta_p_y = 0.01;
                    const double delta_p_z = 0.01;
                    for(double p_x = lower_bound[0]+delta_p_x/2.0; p_x <= upper_bound[0]; p_x += delta_p_x){
                        for(double p_y = lower_bound[1]+delta_p_y/2.0; p_y <= upper_bound[1]; p_y += delta_p_y){
                            for(double p_z = lower_bound[2]+delta_p_z/2.0; p_z <= upper_bound[2]; p_z += delta_p_z){
                                // Make sure both functions produce the same result.
                                PointType test_point = {p_x, p_y, p_z};
                                bool res1 = p_trimmed_domain->IsInsideTrimmedDomain(test_point);
                                bool res2 = brep_operator.IsInside(test_point);
                                QuESo_CHECK_EQUAL(res1, res2);
                            }
                        }
                    }
                    num_of_trimmed_elements++;
                }
            }
        }
    }
    QuESo_CHECK_EQUAL(num_of_trimmed_elements, 701);
}

BOOST_AUTO_TEST_CASE(BunnyPointClassifierOnTrimmedDomainTest) {
    QuESo_INFO << "Testing :: Test Point Classifier On Trimmed Domain:: Bunny Point Classifier" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/stanford_bunny.stl");

    Parameters params( {Component("lower_bound_xyz", PointType(0.0, 0.0, 0.0)),
                        Component("upper_bound_xyz", PointType(1.0, 1.0, 1.0)),
                        Component("lower_bound_uvw", PointType(0.0, 0.0, 0.0)),
                        Component("upper_bound_uvw", PointType(1.0, 1.0, 1.0)),
                        Component("number_of_elements", Vector3i(1, 1, 1)),
                        Component("min_element_volume_ratio", 0.0) });

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh, params);

    const double delta_x = 10;
    const double delta_y = 10;
    const double delta_z = 10;

    IndexType num_of_trimmed_elements = 0;
    for(double x = -24; x <= 85; x += delta_x){
        for(double y = -43; y <= 46; y += delta_y){
            for(double z = 5; z <= 115; z += delta_z){
                Vector3d lower_bound = {x, y, z};
                Vector3d upper_bound = {x+delta_x, y+delta_y, z+delta_z};
                if( brep_operator.GetIntersectionState(lower_bound, upper_bound) == IntersectionStatus::Trimmed){
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound);
                    const double delta_p_x = 2;
                    const double delta_p_y = 2;
                    const double delta_p_z = 2;
                    for(double p_x = lower_bound[0]+delta_p_x/2.0; p_x <= upper_bound[0]; p_x += delta_p_x){
                        for(double p_y = lower_bound[1]+delta_p_y/2.0; p_y <= upper_bound[1]; p_y += delta_p_y){
                            for(double p_z = lower_bound[2]+delta_p_z/2.0; p_z <= upper_bound[2]; p_z += delta_p_z){
                                // Make sure both functions produce the same result.
                                PointType test_point = {p_x, p_y, p_z};
                                bool res1 = p_trimmed_domain->IsInsideTrimmedDomain(test_point);
                                bool res2 = brep_operator.IsInside(test_point);
                                QuESo_CHECK_EQUAL(res1, res2);
                            }
                        }
                    }
                    num_of_trimmed_elements++;
                }
            }
        }
    }
    QuESo_CHECK_EQUAL(num_of_trimmed_elements, 381);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso