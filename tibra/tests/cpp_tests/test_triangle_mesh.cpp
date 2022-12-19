// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "containers/triangle_mesh.h"
#include "io/io_utilities.h"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( TriangleMeshTestSuite )

BOOST_AUTO_TEST_CASE(TriangleMeshIOTest) {
    std::cout << "Testing :: Test Triangle Mesh :: Test IO" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    // Make basic check
    BOOST_CHECK(triangle_mesh.Check());

    BOOST_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 888);
    // Check surface area
    double surface_area = 0.0;

    for( int triangle_id = 0; triangle_id < triangle_mesh.NumOfTriangles(); ++triangle_id){
        surface_area += triangle_mesh.Area(triangle_id);
    }
    BOOST_CHECK_CLOSE(surface_area, 69.11212872984862, 1e-10);
}


BOOST_AUTO_TEST_CASE(TriangleMeshRefineTest) {
    std::cout << "Testing :: Test Triangle Mesh :: Test Refine" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/elephant.stl");

    // Make basic check
    BOOST_CHECK(triangle_mesh.Check());
    double initial_area = 0.0;
    Vector3d initial_weighted_normal(0.0, 0.0, 0.0);
    BOOST_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 5558);
    for( int triangle_id = 0; triangle_id < triangle_mesh.NumOfTriangles(); ++triangle_id){
        initial_area += triangle_mesh.Area(triangle_id);
        auto p_ips = triangle_mesh.GetIPsGlobal(triangle_id, 1);
        for( auto& ip : (*p_ips) ){
            initial_weighted_normal += triangle_mesh.Normal(triangle_id)*ip.GetWeight();
        }
    }


    triangle_mesh.Refine(20000);
    double new_area = 0.0;
    Vector3d new_weighted_normal(0.0, 0.0, 0.0);
    BOOST_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 27911);
    for( int triangle_id = 0; triangle_id < triangle_mesh.NumOfTriangles(); ++triangle_id){
        new_area += triangle_mesh.Area(triangle_id);
        auto p_ips = triangle_mesh.GetIPsGlobal(triangle_id, 1);
        for( auto& ip : (*p_ips) ){
            //std::cout << ip.GetWeight() << std::endl;
            new_weighted_normal += triangle_mesh.Normal(triangle_id)*ip.GetWeight();
        }
    }

    double weighted_normal_error = (initial_weighted_normal- new_weighted_normal).Norm();
    BOOST_CHECK_LT( weighted_normal_error, 1e-14);

    BOOST_CHECK_LT( std::abs(new_area-initial_area)/new_area, 1e-10 );

}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra