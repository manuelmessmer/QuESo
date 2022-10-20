// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <chrono>

#include "embedding/clipper.h"

namespace Testing{

BOOST_AUTO_TEST_SUITE( ClipperTestSuite )

BOOST_AUTO_TEST_CASE(ClipperTest1) {
    std::cout << "Testing :: Test Touching Cube 1" << std::endl;

    // Read mesh from STL file
    // TriangleMesh triangle_mesh{};
    // IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

}

BOOST_AUTO_TEST_SUITE_END()

} // End TouchingCubeTest1
