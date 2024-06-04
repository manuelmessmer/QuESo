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

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/containers/element_container.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/brep_operator.h"
#include "queso/io/io_utilities.h"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( PointClassifierOnTrimmedDomainTestSuite )

BOOST_AUTO_TEST_CASE(PointClassifierOnTrimmedDomainTestSuite) {

    QuESo_INFO << "Testing :: Test Point Classifier On Trimmed Domain:: Cylinder Point Classifier" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cylinder.stl");

    Parameters params( {Component("lower_bound_xyz", PointType{-1.5, -1.5, -1.0}),
                        Component("upper_bound_xyz", PointType{1.5, 1.5, 12}),
                        Component("number_of_elements", Vector3i{6, 6, 26}),
                        Component("min_element_volume_ratio", 0.0) });

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    Mapper mapper(params);
    IndexType num_of_trimmed_elements = 0;
    for( IndexType i = 0; i < mapper.NumberOfElements(); ++i){
        const auto bounding_box = mapper.GetBoundingBoxXYZFromIndex(i);
        const Vector3d lower_bound = bounding_box.first;
        const Vector3d upper_bound = bounding_box.second;
        if( brep_operator.GetIntersectionState(lower_bound, upper_bound) == IntersectionStatus::Trimmed){
            auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound, min_vol_ratio, min_num_triangles);
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
    QuESo_CHECK_EQUAL(num_of_trimmed_elements, 240);
}

BOOST_AUTO_TEST_CASE(CubePointClassifierOnTrimmedDomainTest) {
    QuESo_INFO << "Testing :: Test Point Classifier On Trimmed Domain:: Cube with Cavity Point Classifier" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cube_with_cavity.stl");

    Parameters params( {Component("lower_bound_xyz", PointType{0.0, 0.0, 0.0}),
                        Component("upper_bound_xyz", PointType{1.0, 1.0, 1.0}),
                        Component("number_of_elements", Vector3i{1, 1, 1}),
                        Component("min_element_volume_ratio", 0.0) });

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

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
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound, min_vol_ratio, min_num_triangles);
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

    Parameters params( {Component("lower_bound_xyz", PointType{-0.4, -0.6, -0.35}),
                        Component("upper_bound_xyz", PointType{0.4, 0.6, 0.35}),
                        Component("number_of_elements", Vector3i{16, 24, 14}),
                        Component("min_element_volume_ratio", 0.0) });

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    Mapper mapper(params);
    IndexType num_of_trimmed_elements = 0;
    for( IndexType i = 0; i < mapper.NumberOfElements(); ++i){
        const BoundingBoxType bounding_box = mapper.GetBoundingBoxXYZFromIndex(i);
        Vector3d lower_bound = bounding_box.first;
        Vector3d upper_bound = bounding_box.second;
        if( brep_operator.GetIntersectionState(lower_bound, upper_bound) == IntersectionStatus::Trimmed){
            auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound, min_vol_ratio, min_num_triangles);
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
    QuESo_CHECK_EQUAL(num_of_trimmed_elements, 701);
}

BOOST_AUTO_TEST_CASE(BunnyPointClassifierOnTrimmedDomainTest) {
    QuESo_INFO << "Testing :: Test Point Classifier On Trimmed Domain:: Bunny Point Classifier" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/stanford_bunny.stl");

    Parameters params( {Component("lower_bound_xyz", PointType{-24, -43, 5}),
                        Component("upper_bound_xyz", PointType{85, 46, 115}),
                        Component("lower_bound_uvw", PointType{0.0, 0.0, 0.0}),
                        Component("upper_bound_uvw", PointType{1.0, 1.0, 1.0}),
                        Component("number_of_elements", Vector3i{11, 9, 12})});

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    Mapper mapper(params);
    IndexType num_of_trimmed_elements = 0;
    for( IndexType i = 0; i < mapper.NumberOfElements(); ++i){
        const BoundingBoxType bounding_box = mapper.GetBoundingBoxXYZFromIndex(i);
        Vector3d lower_bound = bounding_box.first;
        Vector3d upper_bound = bounding_box.second;
        if( brep_operator.GetIntersectionState(lower_bound, upper_bound) == IntersectionStatus::Trimmed){
            auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound, min_vol_ratio, min_num_triangles);
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
    QuESo_CHECK_EQUAL(num_of_trimmed_elements, 419);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso