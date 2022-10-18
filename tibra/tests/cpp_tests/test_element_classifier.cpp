// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <chrono>

#include "geometries/triangle_mesh.h"
#include "io/io_utilities.h"
#include "embedding/geometrical_entity_classifier.h"

namespace Testing{

BOOST_AUTO_TEST_SUITE( ElementClassifierTestSuite )

BOOST_AUTO_TEST_CASE(CylinderElementClassifierTest) {
    std::cout << "Testing :: Test AABB Tree :: Ray Intersection" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    GeometricalEntityClassifier classifier(triangle_mesh);

    TriangleMesh::Vector3d lower_bound = {0.0, 0.0, 9.0};
    TriangleMesh::Vector3d upper_bound = {2.0, 2.0, 11.0};
    classifier.GetIntersectionState(lower_bound, upper_bound);


}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing