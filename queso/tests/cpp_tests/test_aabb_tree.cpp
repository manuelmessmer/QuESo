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
//// STL includes
#include <set>
//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/io/io_utilities.h"
#include "queso/embedding/aabb_primitive.h"
#include "queso/embedding/aabb_tree.h"

#include "queso/tests/cpp_tests/global_config.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( AABB_treeTestSuite )

BOOST_AUTO_TEST_CASE(TouchingCubeTest1) {
    QuESo_INFO << "Testing :: Test Touching Cube 1" << std::endl;

    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/cube_with_cavity.stl");

    // Build aabb tree
    AABB_tree tree(triangle_mesh);

    Vector3d lower_bound = {-2.0, -2, -2};
    Vector3d upper_bound = {-1.5, 2, 2};

    AABB_primitive aabb(lower_bound, upper_bound);
    auto results = tree.Query(aabb);

    std::vector<IndexType> intersected_triangles{};

    for( auto r : results){
        const auto& p1 = triangle_mesh.P1(r);
        const auto& p2 = triangle_mesh.P2(r);
        const auto& p3 = triangle_mesh.P3(r);
        // If tolerance>=0 intersection is not detected.
        const double tolerance_1 = 1e-8;
        QuESo_CHECK( !aabb.intersect(p1, p2, p3, tolerance_1) );
        // If tolerance=0 intersection is detected.
        const double tolerance_2 = 0.0;
        QuESo_CHECK( aabb.intersect(p1, p2, p3, tolerance_2) );
    }
    QuESo_CHECK_EQUAL(intersected_triangles.size(), 0);

} // End TouchingCubeTest1


BOOST_AUTO_TEST_CASE(TouchingCubeTest2) {
    QuESo_INFO << "Testing :: Test Touching Cube 2" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/cube_with_cavity.stl");

    // Build aabb tree
    AABB_tree tree(triangle_mesh);

    Vector3d lower_bound = {1.5, -2, -2};
    Vector3d upper_bound = {5.0, 2, 2};

    AABB_primitive aabb(lower_bound, upper_bound);
    auto results = tree.Query(aabb);

    std::vector<IndexType> intersected_triangles{};

    for( auto r : results){
        const auto& p1 = triangle_mesh.P1(r);
        const auto& p2 = triangle_mesh.P2(r);
        const auto& p3 = triangle_mesh.P3(r);
        // If tolerance>=0 intersection is not detected.
        const double tolerance_1 = 1e-8;
        QuESo_CHECK( !aabb.intersect(p1, p2, p3, tolerance_1) );
        // If tolerance=0 intersection is detected.
        const double tolerance_2 = 0.0;
        QuESo_CHECK( aabb.intersect(p1, p2, p3, tolerance_2) );
    }
    QuESo_CHECK_EQUAL(intersected_triangles.size(), 0);

} // End TouchingCubeTest2

BOOST_AUTO_TEST_CASE(TouchingTriangleTest) {
    QuESo_INFO << "Testing :: Test Touching Triangle " << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};

    // Construct inclined triangle.
    triangle_mesh.AddVertex({0.0,0.0,0.0});
    triangle_mesh.AddVertex({-1.0,1.0,0.0});
    triangle_mesh.AddVertex({-1.0,1.0,1.0});
    triangle_mesh.AddTriangle({0, 1, 2});

    // Build aabb tree
    AABB_tree tree(triangle_mesh);

    // Lower bound touched triangle.
    Vector3d lower_bound = {-0.5, 0.5, 0.1};
    Vector3d upper_bound = {1.0, 1.0, 2.0};

    AABB_primitive aabb(lower_bound, upper_bound);

    // This test is conservative and must detect triangle.
    auto results = tree.Query(aabb);

    std::vector<IndexType> intersected_triangles{};

    for( auto r : results){
        const auto& p1 = triangle_mesh.P1(r);
        const auto& p2 = triangle_mesh.P2(r);
        const auto& p3 = triangle_mesh.P3(r);

        // If tolerance>=0 intersection is not detected.
        const double tolerance_1 = 1e-8;
        QuESo_CHECK( !aabb.intersect(p1, p2, p3, tolerance_1) );

        // If tolerance=0 intersection is detected.
        const double tolerance_2 = 0.0;
        QuESo_CHECK( aabb.intersect(p1, p2, p3, tolerance_2) );
    }
} // End TouchingTriangleTest

BOOST_AUTO_TEST_CASE(CylinderFindIntersectedTrianglesTest) {
    QuESo_INFO << "Testing :: Test Find Intersected Triangles :: Cylinder" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/cylinder.stl");

    // Build aabb tree
    AABB_tree tree(triangle_mesh);

    // Open results file. Results are checked with CGAL
    std::ifstream myfile(base_dir + "/results/aabb_cylinder.txt");

    IndexType number_trimmed_elements = 0;
    for(double xx = -1.5; xx <= 1.5; xx += 0.1){
        for(double yy = -1.5; yy <= 1.5; yy += 0.1){
            for(double zz = -1; zz <= 11; zz += 0.1){

                Vector3d lower_bound = {xx, yy, zz};
                Vector3d upper_bound = {xx+0.1, yy+0.1, zz+0.1};

                AABB_primitive aabb(lower_bound, upper_bound);
                auto results = tree.Query(aabb);

                TriangleMesh new_mesh{};

                std::set<int> indices{};
                std::map<int,int> index_map{};

                std::vector<IndexType> intersected_triangles{};

                for( auto r : results){
                    const auto& p1 = triangle_mesh.P1(r);
                    const auto& p2 = triangle_mesh.P2(r);
                    const auto& p3 = triangle_mesh.P3(r);
                    const double tolerance = 0.0;
                    if( aabb.intersect(p1, p2, p3, tolerance) ){
                        intersected_triangles.push_back(r);
                    }
                }

                if( intersected_triangles.size() > 0 ){
                    std::string line;
                    // Read header: Element number
                    getline(myfile,line);
                    for( auto& rr : intersected_triangles ){
                        getline (myfile,line);
                        const int triangle_id_ref = std::stoi(line);
                        QuESo_CHECK_EQUAL(triangle_id_ref, static_cast<int>(rr));
                    }
                    number_trimmed_elements++;
                }
            }
        }
    }
    QuESo_CHECK_EQUAL(number_trimmed_elements, 8212);

    myfile.close();
} // End CylinderFindIntersectedTrianglesTest

BOOST_AUTO_TEST_CASE(ElephantFindIntersectedTrianglesTest) {
    QuESo_INFO << "Testing :: Test Find Intersected Triangles :: Elephant" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/elephant.stl");

    // Build aabb tree.
    AABB_tree tree(triangle_mesh);

    // Open results file. Results are verified with CGAL. Touching triangles are not considered as intersected.
    std::ifstream myfile(base_dir + "/results/aabb_elephant.txt");

    IndexType number_trimmed_elements = 0;
    for(double xx = -0.4; xx <= 0.4; xx += 0.02){
        for(double yy = -0.6; yy <= 0.6; yy += 0.02){
            for(double zz = -0.35; zz <= 0.35; zz += 0.02){
                Vector3d lower_bound = {xx, yy, zz};
                Vector3d upper_bound = {xx+0.02, yy+0.02, zz+0.02};

                AABB_primitive aabb(lower_bound, upper_bound);
                auto results = tree.Query(aabb);

                std::vector<IndexType> intersected_triangles{};
                for( auto r : results){
                    const auto& p1 = triangle_mesh.P1(r);
                    const auto& p2 = triangle_mesh.P2(r);
                    const auto& p3 = triangle_mesh.P3(r);
                    const double tolerance = 0.0;
                    if( aabb.intersect(p1, p2, p3, tolerance) ){
                        intersected_triangles.push_back(r);
                    }
                }

                if( intersected_triangles.size() > 0 ){
                    std::string line;
                    // Read header: Element number
                    getline(myfile,line);
                    for( auto& rr : intersected_triangles ){
                        getline (myfile,line);
                        const int triangle_id_ref = std::stoi(line);
                        QuESo_CHECK_EQUAL(triangle_id_ref, static_cast<int>(rr));
                    }
                    number_trimmed_elements++;
                }
            }
        }
    }

    QuESo_CHECK_EQUAL(number_trimmed_elements, 4641);
    myfile.close();

} // End ElephantFindIntersectedTrianglesTest

BOOST_AUTO_TEST_CASE(BunnyFindIntersectedTrianglesTest) {
    QuESo_INFO << "Testing :: Test Find Intersected Triangles :: Bunny" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/stanford_bunny.stl");

    // Build aabb tree.
    AABB_tree tree(triangle_mesh);

    // Read results file. Results are verified with CGAL:
    std::ifstream myfile(base_dir + "/results/aabb_bunny.txt");

    IndexType number_trimmed_elements = 0;
    for(double xx = -24; xx <= 85; xx += 2){
        for(double yy = -43; yy <= 46; yy += 2){
            for(double zz = 5; zz <= 115; zz += 2){
                Vector3d lower_bound = {xx, yy, zz};
                Vector3d upper_bound = {xx+2, yy+2, zz+2};

                AABB_primitive aabb(lower_bound, upper_bound);
                auto results = tree.Query(aabb);

                std::vector<IndexType> intersected_triangles{};
                for( auto r : results){
                    const auto& p1 = triangle_mesh.P1(r);
                    const auto& p2 = triangle_mesh.P2(r);
                    const auto& p3 = triangle_mesh.P3(r);
                    const double tolerance = 0.0;
                    if( aabb.intersect(p1, p2, p3, tolerance) ){
                        intersected_triangles.push_back(r);
                    }

                }

                if( intersected_triangles.size() > 0 ){
                    std::string line;
                    // Read header: Element number
                    getline(myfile,line);
                    for( auto& rr : intersected_triangles ){
                        getline (myfile,line);
                        const int triangle_id_ref = std::stoi(line);
                        QuESo_CHECK_EQUAL(triangle_id_ref, static_cast<int>(rr));
                    }
                    number_trimmed_elements++;
                }
            }
        }
    }
    QuESo_CHECK_EQUAL(number_trimmed_elements, 10531);
    myfile.close();

} // End BunnyFindIntersectedTrianglesTest

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso
