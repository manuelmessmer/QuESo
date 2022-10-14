// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <chrono>

#include "geometries/triangle_mesh.h"
#include "io/io_utilities.h"
#include "utilities/aabb_tree.h"

namespace Testing{

BOOST_AUTO_TEST_SUITE( AABB_treeTestSuite )

BOOST_AUTO_TEST_CASE(CylinderTriangleIntersectionTest) {
    std::cout << "Testing :: Test AABB Tree :: Ray Intersection" << std::endl;

    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    AABB_tree tree(triangle_mesh);

    TriangleMesh::Vector3d lower_bound = {0, 0, 4.0};
    TriangleMesh::Vector3d upper_bound = {2, 2, 5.0};

    AABB_primitive aabb(lower_bound, upper_bound);
    auto results = tree.Query(aabb);

    TriangleMesh new_mesh{};

    std::set<int> indices{};
    std::map<int,int> index_map{};

    std::vector<IndexType> intersected_triangles{};
    int count = 0;
    for( auto r : results){
        const auto& p1 = triangle_mesh.P1(r);
        const auto& p2 = triangle_mesh.P2(r);
        const auto& p3 = triangle_mesh.P3(r);
        double t, u, v;

        if( aabb.intersect(p1, p2, p3, t, u, v) ){
            intersected_triangles.push_back(r);
        }
    }

    new_mesh.Copy(intersected_triangles, triangle_mesh);
    IO::WriteMeshToSTL(new_mesh,"test.stl", true);

}



BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing