// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "containers/triangle_mesh.h"
#include "utilities/mesh_utilities.h"
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

std::pair<double,Vector3d> ComputeAreaAndWeightedNormal(const TriangleMesh& rTriangleMesh){
    double area = 0.0;
    Vector3d weighted_normal(0.0, 0.0, 0.0);

    for( int triangle_id = 0; triangle_id < rTriangleMesh.NumOfTriangles(); ++triangle_id){
        area += rTriangleMesh.Area(triangle_id);
        auto p_ips = rTriangleMesh.pGetIPsGlobal(triangle_id, 1);
        for( auto& ip : (*p_ips) ){
            weighted_normal += rTriangleMesh.Normal(triangle_id)*ip.GetWeight();
        }
    }
    return std::make_pair(area, weighted_normal);
}

BOOST_AUTO_TEST_CASE(TriangleMeshRefineTest) {
    std::cout << "Testing :: Test Triangle Mesh :: Test Refine" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/elephant.stl");

    // Make basic check
    BOOST_CHECK(triangle_mesh.Check());
    BOOST_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 5558);
    auto init_values = ComputeAreaAndWeightedNormal(triangle_mesh);

    MeshUtilities::Refine(triangle_mesh, 20000);
    BOOST_CHECK(triangle_mesh.Check());
    BOOST_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 27911);
    auto new_values = ComputeAreaAndWeightedNormal(triangle_mesh);

    double weighted_normal_error = (init_values.second - new_values.second).Norm();
    BOOST_CHECK_LT( weighted_normal_error, 1e-14);
    BOOST_CHECK_LT( std::abs(init_values.first - new_values.first)/init_values.first, 1e-10 );
}

BOOST_AUTO_TEST_CASE(TriangleMeshAppendTest) {
    std::cout << "Testing :: Test Triangle Mesh :: Test Append" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/stanford_bunny.stl");

    // Make basic check
    BOOST_CHECK(triangle_mesh.Check());
    BOOST_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 112402);
    auto init_values = ComputeAreaAndWeightedNormal(triangle_mesh);

    TriangleMesh new_mesh{};
    MeshUtilities::Append(new_mesh, triangle_mesh);
    BOOST_CHECK(new_mesh.Check());
    BOOST_CHECK_EQUAL(new_mesh.NumOfTriangles(), 112402);
    auto new_values = ComputeAreaAndWeightedNormal(new_mesh);

    double weighted_normal_error = (init_values.second - new_values.second).Norm();
    BOOST_CHECK_LT( weighted_normal_error, 1e-14);
    BOOST_CHECK_LT( std::abs(init_values.first - new_values.first)/init_values.first, 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra