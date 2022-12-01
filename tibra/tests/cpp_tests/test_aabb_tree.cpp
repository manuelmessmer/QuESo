// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <set>
#include <chrono>

#include "containers/triangle_mesh.h"
#include "io/io_utilities.h"
#include "embedding/aabb_tree.h"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( AABB_treeTestSuite )

BOOST_AUTO_TEST_CASE(TouchingCubeTest1) {
    std::cout << "Testing :: Test Touching Cube 1" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

    // Build aabb tree
    AABB_tree tree(triangle_mesh);

    Vector3d lower_bound = {-2.0, -2, -2};
    Vector3d upper_bound = {-1.5, 2, 2};

    AABB_primitive aabb(lower_bound, upper_bound);
    auto results = tree.Query(aabb);

    std::vector<IndexType> intersected_triangles{};
    int count = 0;

    for( auto r : results){
        const auto& p1 = triangle_mesh.P1(r);
        const auto& p2 = triangle_mesh.P2(r);
        const auto& p3 = triangle_mesh.P3(r);
        // If tolerance>=0 intersection is not detected.
        const double tolerance_1 = 1e-8;
        BOOST_CHECK( !aabb.intersect(p1, p2, p3, tolerance_1) );
        // If tolerance=0 intersection is detected.
        const double tolerance_2 = 0.0;
        BOOST_CHECK( aabb.intersect(p1, p2, p3, tolerance_2) );
    }
    BOOST_CHECK_EQUAL(intersected_triangles.size(), 0);

} // End TouchingCubeTest1


BOOST_AUTO_TEST_CASE(TouchingCubeTest2) {
    std::cout << "Testing :: Test Touching Cube 2" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

    // Build aabb tree
    AABB_tree tree(triangle_mesh);

    Vector3d lower_bound = {1.5, -2, -2};
    Vector3d upper_bound = {5.0, 2, 2};

    AABB_primitive aabb(lower_bound, upper_bound);
    auto results = tree.Query(aabb);

    std::vector<IndexType> intersected_triangles{};
    int count = 0;

    for( auto r : results){
        const auto& p1 = triangle_mesh.P1(r);
        const auto& p2 = triangle_mesh.P2(r);
        const auto& p3 = triangle_mesh.P3(r);
        // If tolerance>=0 intersection is not detected.
        const double tolerance_1 = 1e-8;
        BOOST_CHECK( !aabb.intersect(p1, p2, p3, tolerance_1) );
        // If tolerance=0 intersection is detected.
        const double tolerance_2 = 0.0;
        BOOST_CHECK( aabb.intersect(p1, p2, p3, tolerance_2) );
    }
    BOOST_CHECK_EQUAL(intersected_triangles.size(), 0);

} // End TouchingCubeTest2

BOOST_AUTO_TEST_CASE(TouchingTriangleTest) {
    std::cout << "Testing :: Test Touching Triangle " << std::endl;

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
    int count = 0;

    for( auto r : results){
        const auto& p1 = triangle_mesh.P1(r);
        const auto& p2 = triangle_mesh.P2(r);
        const auto& p3 = triangle_mesh.P3(r);

        // If tolerance>=0 intersection is not detected.
        const double tolerance_1 = 1e-8;
        BOOST_CHECK( !aabb.intersect(p1, p2, p3, tolerance_1) );

        // If tolerance=0 intersection is detected.
        const double tolerance_2 = 0.0;
        BOOST_CHECK( aabb.intersect(p1, p2, p3, tolerance_2) );
    }
} // End TouchingTriangleTest

BOOST_AUTO_TEST_CASE(CylinderFindIntersectedTrianglesTest) {
    std::cout << "Testing :: Find Intersected Triangles :: Cylinder" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    // Build aabb tree
    AABB_tree tree(triangle_mesh);

    // Open results file. Results are checked with CGAL
    std::ifstream myfile("tibra/tests/cpp_tests/results/aabb_cylinder.txt");

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
                int count = 0;

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
                        BOOST_CHECK_EQUAL(triangle_id_ref, rr);
                    }
                    number_trimmed_elements++;
                }
            }
        }
    }
    BOOST_CHECK_EQUAL(number_trimmed_elements, 8212);

    myfile.close();
} // End CylinderFindIntersectedTrianglesTest

BOOST_AUTO_TEST_CASE(ElephantFindIntersectedTrianglesTest) {
    std::cout << "Testing :: Find Intersected Triangles :: Elephant" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/elephant.stl");

    // Build aabb tree.
    AABB_tree tree(triangle_mesh);

    // Open results file. Results are verified with CGAL. Touching triangles are not considered as intersected.
    std::ifstream myfile("tibra/tests/cpp_tests/results/aabb_elephant.txt");

    IndexType number_trimmed_elements = 0;
    for(double xx = -0.4; xx <= 0.4; xx += 0.02){
        for(double yy = -0.6; yy <= 0.6; yy += 0.02){
            for(double zz = -0.35; zz <= 0.35; zz += 0.02){
                Vector3d lower_bound = {xx, yy, zz};
                Vector3d upper_bound = {xx+0.02, yy+0.02, zz+0.02};

                AABB_primitive aabb(lower_bound, upper_bound);
                auto results = tree.Query(aabb);

                std::vector<IndexType> intersected_triangles{};
                int count = 0;
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
                        BOOST_CHECK_EQUAL(triangle_id_ref, rr);
                    }
                    number_trimmed_elements++;
                }
            }
        }
    }

    BOOST_CHECK_EQUAL(number_trimmed_elements, 4641);
    myfile.close();

} // End ElephantFindIntersectedTrianglesTest

BOOST_AUTO_TEST_CASE(BunnyFindIntersectedTrianglesTest) {
    std::cout << "Testing :: Find Intersected Triangles :: Bunny" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/stanford_bunny.stl");

    // Build aabb tree.
    AABB_tree tree(triangle_mesh);

    // Read results file. Results are verified with CGAL:
    std::ifstream myfile("tibra/tests/cpp_tests/results/aabb_bunny.txt");

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
                        BOOST_CHECK_EQUAL(triangle_id_ref, rr);
                    }
                    number_trimmed_elements++;
                }
            }
        }
    }
    BOOST_CHECK_EQUAL(number_trimmed_elements, 10531);
    myfile.close();

} // End BunnyFindIntersectedTrianglesTest

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra


// #include "containers/element.h"
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Mesh_polyhedron_3.h>
// #include <CGAL/Surface_mesh.h>
// #include <CGAL/boost/graph/IO/STL.h>
// #include <CGAL/AABB_face_graph_triangle_primitive.h>
// #include <CGAL/Polygon_mesh_processing/locate.h>
// #include <CGAL/AABB_tree.h>
// #include <CGAL/Polygon_mesh_processing/intersection.h>
// #include <CGAL/Polygon_mesh_processing/corefinement.h>
// #include <CGAL/Polygon_mesh_processing/clip.h>

// BOOST_AUTO_TEST_CASE(ElephantTriangleIntersectionTest) {
//     std::cout << "Testing :: Test AABB Tree :: Ray Intersection" << std::endl;


//     //CGAl test
//     typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//     typedef K::Point_3 Point_3;
//     typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
//     typedef CGAL::AABB_face_graph_triangle_primitive<SurfaceMeshType>           AABBFaceGraphsPrimitivesType;
//     typedef CGAL::AABB_traits<K, AABBFaceGraphsPrimitivesType>                  AABBFaceGraphsTraits;
//     typedef CGAL::AABB_tree<AABBFaceGraphsTraits>                               AABBTreeType;
//     typedef Element::PositionType PositionType;

//     SurfaceMeshType mPolyhedron;
//     CGAL::IO::read_STL("tibra/tests/cpp_tests/data/elephant.stl", mPolyhedron);
//     AABBTreeType tree_cgal{};
//     CGAL::Polygon_mesh_processing::build_AABB_tree(mPolyhedron, tree_cgal);

//     TriangleMesh triangle_mesh{};
//     // Read mesh from STL file
//     IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/elephant.stl");

//     AABB_tree tree(triangle_mesh);

//     std::ifstream myfile;
//     myfile.open ("tibra/tests/cpp_tests/results/aabb_elephant.txt");
//     int count_el = 0;
//     for(double xx = -0.4; xx <= 0.4; xx += 0.02){
//         for(double yy = -0.6; yy <= 0.6; yy += 0.02){
//             for(double zz = -0.35; zz <= 0.35; zz += 0.02){

//             Vector3d lower_bound = {xx, yy, zz};
//             Vector3d upper_bound = {xx+0.02, yy+0.02, zz+0.02};

//             const Point_3 point1(lower_bound[0], lower_bound[1], lower_bound[2]);
//             const Point_3 point2(upper_bound[0], upper_bound[1], upper_bound[2]);
//             const CGAL::Iso_cuboid_3<K> tmp_cuboid( point1, point2, 0);

//             typedef boost::optional<AABBTreeType::Intersection_and_primitive_id<CGAL::Iso_cuboid_3<K>>::Type> Ray_intersection;

//             std::vector<AABBTreeType::Primitive_id> primitives;

//             auto p = tree_cgal.all_intersected_primitives(tmp_cuboid, std::back_inserter(primitives));

//             PositionType positions = mPolyhedron.points();

//             std::vector<IndexType> actual_cgal_intersections{};

//             for( auto pri : primitives){
//                  Point_3 p1{};
//                 Point_3 p2{};
//                 Point_3 p3{};
//                 std::array<Point_3,3> points_ = {p1, p2, p3};
//                 int ii = 0;
//                 for( auto vh : vertices_around_face(halfedge(pri, mPolyhedron), mPolyhedron) ){
//                     points_[ii] = positions[vh];
//                     ii++;
//                 }
//                 // //const K::Triangle_3 tri(points_[0], points_[1], points_[2]);
//                 p1 = points_[0];
//                 p2 = points_[1];
//                 p3 = points_[2];
//                 // const auto& p1 = triangle_mesh.P1(pri.idx());
//                 // const auto& p2 = triangle_mesh.P2(pri.idx());
//                 // const auto& p3 = triangle_mesh.P3(pri.idx());

//                 double x_min = std::min({p1[0], p2[0], p3[0] });
//                 double y_min = std::min({p1[1], p2[1], p3[1] });
//                 double z_min = std::min({p1[2], p2[2], p3[2] });

//                 double x_max = std::max({p1[0], p2[0], p3[0] });
//                 double y_max = std::max({p1[1], p2[1], p3[1] });
//                 double z_max = std::max({p1[2], p2[2], p3[2] });

//                 if(    x_max <= lower_bound[0] || x_min >= upper_bound[0]
//                     || y_max <= lower_bound[1] || y_min >= upper_bound[1]
//                     || z_max <= lower_bound[2] || z_min >= upper_bound[2] )
//                 {
//                     std::cout << std::setprecision(15) << std::endl;
//                     std::cout << "xx: " << lower_bound[0] << std::endl;
//                     std::cout << "yy: " << lower_bound[1] << std::endl;
//                     std::cout << "zz: " << lower_bound[2] << std::endl;
//                     std::cout << "upper_bound 0: " << upper_bound[0] << std::endl;
//                     std::cout << "upper_bound 1: " << upper_bound[1] << std::endl;
//                     std::cout << "upper_bound 2: " << upper_bound[2] << std::endl;
//                     // std::cout << p1 << std::endl;
//                     // std::cout << p2 << std::endl;
//                     // std::cout << p3 << std::endl;
//                 }
//                 else {
//                 }

//                 actual_cgal_intersections.push_back( pri.idx() );

//             }

//             AABB_primitive aabb(lower_bound, upper_bound);
//             auto results = tree.Query(aabb);

//             TriangleMesh new_mesh{};

//             std::set<int> indices{};
//             std::map<int,int> index_map{};

//             std::vector<IndexType> intersected_triangles{};
//             int count = 0;
//             for( auto r : results){
//                 const auto& p1 = triangle_mesh.P1(r);
//                 const auto& p2 = triangle_mesh.P2(r);
//                 const auto& p3 = triangle_mesh.P3(r);
//                 double t, u, v;

//                     //actual_cgal_intersections.push_back( pri.idx() );
//                 if( aabb.intersect(p1, p2, p3, t, u, v) ){
//                     intersected_triangles.push_back(r);
//                 }

//             }

//             if( intersected_triangles.size() > 0 ){
//                 std::string line;
//                 // Read header: Element number
//                 getline(myfile,line);
//                 for( auto& rr : intersected_triangles ){
//                     getline (myfile,line);
//                     const int triangle_id_ref = std::stoi(line);
//                     BOOST_CHECK_EQUAL(triangle_id_ref, rr);
//                 }
//             }

//             // if( intersected_triangles.size() > 0 )
//             //     std::cout << "size: " << intersected_triangles.size() << std::endl;
//             int cc = 0;
//             std::sort(intersected_triangles.begin(), intersected_triangles.end());
//             std::sort(actual_cgal_intersections.begin(), actual_cgal_intersections.end() );
//             //std::cout << intersected_triangles[0] << std::endl;
//             if( !(actual_cgal_intersections.size() == intersected_triangles.size())){
//                 std::cout << "x: " << xx << ", " << yy << ", " << zz << std::endl;
//                 std::cout << "x: " << xx+0.1 << ", " << yy+0.1 << ", " << zz+0.1 << std::endl;
//             }
//             //BOOST_CHECK_EQUAL(actual_cgal_intersections.size(), intersected_triangles.size() );
//             if( intersected_triangles.size() > 0){
//                 std::cout << actual_cgal_intersections.size() << ": " << intersected_triangles.size() << std::endl;
//             }
//             for( int i = 0; i < actual_cgal_intersections.size(); ++i){
//                 //std::cout << actual_cgal_intersections[i] << ", " << intersected_triangles[i] << std::endl;
//                 if( !(actual_cgal_intersections.size() == intersected_triangles.size())){
//                     std::cout << "x: " << xx << ", " << yy << ", " << zz << std::endl;
//                     std::cout << "x: " << xx+0.1 << ", " << yy+0.1 << ", " << zz+0.1 << std::endl;
//                     const auto p1 = triangle_mesh.P1(actual_cgal_intersections[i]);
//                     std::cout << p1[0] << ", " << p1[1] << ", " << p1[2] << ", " << std::endl;
//                     const auto p2 = triangle_mesh.P2(actual_cgal_intersections[i]);
//                     std::cout << p2[0] << ", " << p2[1] << ", " << p2[2] << ", " << std::endl;
//                     const auto p3 = triangle_mesh.P3(actual_cgal_intersections[i]);
//                     std::cout << p3[0] << ", " << p3[1] << ", " << p3[2] << ", " << std::endl;
//                 }
//                 BOOST_CHECK_EQUAL(actual_cgal_intersections[i], intersected_triangles[i]);
//             }

//             }
//         }
//     }
//     // std::cout << intersected_triangles[cc] << std::endl;
//     // new_mesh.Copy(intersected_triangles, triangle_mesh);
//     // IO::WriteMeshToSTL(new_mesh,"test.stl", true);
//     myfile.close();

// }