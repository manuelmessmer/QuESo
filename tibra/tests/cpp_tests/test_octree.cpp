// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// STL includes

//// Project includes
#include "containers/triangle_mesh.hpp"
#include "io/io_utilities.h"
#include "embedding/brep_operator.h"
#include "embedding/octree.h"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( Octree_TestSuite )

BOOST_AUTO_TEST_CASE(TouchingCubeTest1) {
    TIBRA_INFO << "Testing :: Test Touching Cube 1" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

    Parameters params{};
    // Build aabb tree
    //AABB_tree tree(triangle_mesh);
    BRepOperator brep_operator(triangle_mesh, params);
    auto p_trimmed_domain = brep_operator.GetTrimmedDomain({-2.0, -2, -2},{-1.3, -1.3, -1.3});
    auto mesh = p_trimmed_domain->GetTriangleMesh();

    Octree octree(p_trimmed_domain.get());
    octree.Refine(4, 6);

    std::cout << "get: " << octree.Volume() << std::endl;

    IO::WriteMeshToSTL(mesh, "mesh.stl", true);
    // Vector3d lower_bound = {-2.0, -2, -2};
    // Vector3d upper_bound = {-1.5, 2, 2};


} // End TouchingCubeTest1


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra

