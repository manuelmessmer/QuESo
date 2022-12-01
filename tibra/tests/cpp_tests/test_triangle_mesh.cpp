// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "containers/triangle_mesh.h"
#include "io/io_utilities.h"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( TriangleMeshTestSuite )

BOOST_AUTO_TEST_CASE(TriangleMeshTest1) {

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

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra