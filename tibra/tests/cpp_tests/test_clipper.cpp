// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <chrono>

#include "geometries/triangle_mesh.h"
#include "io/io_utilities.h"
#include "embedding/brep_operator.h"

namespace Testing{

BOOST_AUTO_TEST_SUITE( ClipperTestSuite )

BOOST_AUTO_TEST_CASE(ClipperTest1) {
    std::cout << "Testing :: Test Touching Cube 1" << std::endl;

    //Read mesh from STL file
    TriangleMesh triangle_mesh{};

    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    TriangleMesh::Vector3d lower_bound = {0.0, 0.0, 0.0};
    TriangleMesh::Vector3d upper_bound = {2, 2, 1};

    auto p_clipped_mesh = brep_operator.ClipTriangleMesh(lower_bound, upper_bound);

    p_clipped_mesh->Check();

    double area = 0.0;
    for( int i = 0; i < p_clipped_mesh->NumOfTriangles(); ++i){
        area += p_clipped_mesh->Area(i);
    }

    // Check surface area of clipped mesh.
    BOOST_CHECK_CLOSE( area, 1.570796327, 1e-2);

}

BOOST_AUTO_TEST_SUITE_END()

} // End TouchingCubeTest1
