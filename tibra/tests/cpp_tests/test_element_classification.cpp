// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <chrono>

#include "geometries/triangle_mesh.h"
#include "io/io_utilities.h"
#include "utilities/element_classification.h"


/// Remove this again
/// CGAL includes
// Domain
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/STL.h>
#include "utilities/intersection_test.h"

namespace Testing{

BOOST_AUTO_TEST_SUITE( ElementClassificationTestSuite )

BOOST_AUTO_TEST_CASE(InsideOutsideTest1) {

    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    std::vector<PointType> rPoints{};


    for(double x = -1.5; x <= 1.5; x += 0.1){
        for(double y = -1.5; y <= 1.5; y += 0.1){
            for(double z = -1; z <= 12; z += 0.1){

                if(x < -1.0
                || x > 1.0
                || y < -1.0
                || y > 1.0
                || z < 0.0
                || z > 10.0)
                {
                //return CGAL::ON_UNBOUNDED_SIDE;
                }
                else {
                    rPoints.push_back( {x, y, z} );

                }

            }
        }
    }



    std::cout << "Size: " << rPoints.size() << '\n';
    std::cout << "start test: \n";
    auto start_time = std::chrono::high_resolution_clock::now();
    auto result = ElementClassification::PointsAreInside(triangle_mesh, rPoints);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "end test: \n";
    for( int i = 0; i < result->size(); ++i){

        double radius = std::sqrt( rPoints[i][0]*rPoints[i][0] + rPoints[i][1]*rPoints[i][1] );
        if( radius < 1.0-1e-14 && rPoints[i][2] >= 0.0 && rPoints[i][2] <= 10.0){
            if( !(*result)[i] )
                std::cout << "Radius: " << radius << ", " << rPoints[i][2] << std::endl;
            BOOST_CHECK((*result)[i]);
        }
        else {
            BOOST_CHECK(!(*result)[i]);
        }
        //std::cout << result[i] << std::endl;
    }

    std::chrono::duration<double> elapsed_time = end_time - start_time;
    std::cout << "TIBRA :: Elapsed Time: " << elapsed_time.count() << std::endl;
    IO::WriteMeshToSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder2.stl", true);

}




BOOST_AUTO_TEST_CASE(InsideOutsideTest2) {

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3;
    typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;

    SurfaceMeshType triangle_mesh{};
    // Read mesh from STL file
    //IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");
    CGAL::IO::read_STL("tibra/tests/cpp_tests/data/cylinder.stl", triangle_mesh);
    std::vector<PointType> rPoints{};

    std::array<double,3> point_a{0.0, 0.0 ,0.0};
    std::array<double,3> point_b{1.0, 1.0 ,1.0};
    IntersectionTest intersection_test(triangle_mesh, point_a, point_b);


    for(double x = -1.5; x <= 1.5; x += 0.1){
        for(double y = -1.5; y <= 1.5; y += 0.1){
            for(double z = -1; z <= 12; z += 0.1){
                rPoints.push_back( {x, y, z} );
            }
        }
    }

    std::cout << "Size: " << rPoints.size() << '\n';
    std::cout << "start test: \n";
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<bool> result(rPoints.size(), false);
    int count = 0;
    for( auto& point : rPoints){
        if( intersection_test.IsInside(point) ){
            result[count] = true;
        }

        count++;
    }
    std::cout << "count " << count << std::endl;
    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "end test: \n";
    for( int i = 0; i < result.size(); ++i){

        double radius = std::sqrt( rPoints[i][0]*rPoints[i][0] + rPoints[i][1]*rPoints[i][1] );
        if( radius < 1.0-1e-14 && rPoints[i][2] >= 0.0 && rPoints[i][2] <= 10.0){
            if( !(result)[i] )
                std::cout << "Radius: " << radius << ", " << rPoints[i][2] << std::endl;
            BOOST_CHECK((result)[i]);
        }
        else {
            BOOST_CHECK(!(result)[i]);
        }
        //std::cout << result[i] << std::endl;
    }

    std::chrono::duration<double> elapsed_time = end_time - start_time;
    std::cout << "TIBRA :: Elapsed Time: " << elapsed_time.count() << std::endl;
    //IO::WriteMeshToSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder2.stl", true);

}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing