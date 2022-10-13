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

BOOST_AUTO_TEST_CASE(RayIntersectionTest) {
    std::cout << "Testing :: Test AABB Tree :: Ray Intersection" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    AABB_tree aabb_tree(triangle_mesh);
    Ray_AABB_primitive ray({0.0, 0.0, 1.0}, {1.0, 1.0, 0.0});
    auto result = aabb_tree.Query(ray);

    std::cout << "Size: " << result.size() << std::endl;
    for( auto r : result ){
        const auto& p1 = triangle_mesh.P1(r);
        const auto& p2 = triangle_mesh.P2(r);
        const auto& p3 = triangle_mesh.P3(r);

        double t, u, v;
        if( ray.intersect(p1, p2, p3, t, u, v) ) {
            std::cout << "id: " << r << std::endl;
            std::array<double,3> pp = {u, v, 0.0};
            auto point_global = triangle_mesh.GetPointGlobalSpace(r, pp);
            std::cout << "t: " << t << " u: " << u  << " v: " << v << std::endl;
            const double sum = u+v;
            if( u < 0.0+1e-15 || v < 0.0+1e-15 || sum > 1.0-1e-15 ){
                std::cout << "on boundary " << std::endl;
            }
            std::cout << "x: " << point_global[0] << " y: " << point_global[1]  << " z: " << point_global[2] << std::endl;
        }
    }
}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing