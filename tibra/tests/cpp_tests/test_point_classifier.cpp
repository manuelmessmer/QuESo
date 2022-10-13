// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <chrono>

#include "geometries/triangle_mesh.h"
#include "io/io_utilities.h"
#include "utilities/geometrical_entity_classifier.h"
#include "geometries/element_container.h"

namespace Testing{

BOOST_AUTO_TEST_SUITE( PointClassifierTestSuite )

BOOST_AUTO_TEST_CASE(CylinderPointClassifierTest) {
    std::cout << "Testing :: Test Geometrical Entity Classifier :: Cylinder Point Classifier" << std::endl;

    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    std::vector<PointType> rPoints{};
    rPoints.reserve(4340000);
    for(double x = -1.5; x <= 1.5; x += 0.03){
        for(double y = -1.5; y <= 1.5; y += 0.03){
            for(double z = -1; z <= 12; z += 0.03){
                rPoints.push_back( {x, y, z} );
            }
        }
    }

    GeometricalEntityClassifier classifier(triangle_mesh);
    //auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<bool> result(rPoints.size(), false);
    int count = 0;
    for( auto& point : rPoints){
        if( classifier.IsInside(point) ){
            result[count] = true;
        }
        count++;
    }
    //auto end_time = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> elapsed_time = end_time - start_time;
    //std::cout << "TIBRA :: Elapsed Time: " << elapsed_time.count() << std::endl;

    for( int i = 0; i < result.size(); ++i){
        double radius = std::sqrt( rPoints[i][0]*rPoints[i][0] + rPoints[i][1]*rPoints[i][1] );
        if( radius < 1.0 && rPoints[i][2] > 0.0 && rPoints[i][2] < 10.0){
            BOOST_CHECK((result)[i]);
        }
        else {
            BOOST_CHECK(!(result)[i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(ElephantPointClassifierTest) {
    std::cout << "Testing :: Test Geometrical Entity Classifier :: Elphant Point Classifier" << std::endl;

    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/elephant.stl");

    GeometricalEntityClassifier classifier(triangle_mesh);

    std::vector<PointType> rPoints{};
    rPoints.reserve(84000);
    for(double x = -0.4; x <= 0.4; x += 0.02){
        for(double y = -0.6; y <= 0.6; y += 0.02){
            for(double z = -0.35; z <= 0.35; z += 0.02){
                rPoints.push_back( {x, y, z} );
            }
        }
    }

    // auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<bool> result(rPoints.size(), false);
    int count = 0;
    for( auto& point : rPoints){
        if( classifier.IsInside(point) ){
            result[count] = true;
        }
        count++;
    }
    // auto end_time = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed_time = end_time - start_time;

    std::vector<bool> result_ref{};
    // Read reference results from file
    std::string line;
    std::ifstream myfile ("tibra/tests/cpp_tests/results/inside_outside_elephant.txt");
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
        result_ref.push_back( std::stoi(line) );
        }
        myfile.close();
    }

    // std::ofstream myfile;
    // myfile.open ("test.txt");
    // Compare results
    for( int i = 0; i < result.size(); ++i){
        BOOST_CHECK_EQUAL(result[i], result_ref[i]);
    }
    // myfile.close();
}


BOOST_AUTO_TEST_CASE(BunnyPointClassifierTest) {
    std::cout << "Testing :: Test Geometrical Entity Classifier :: Bunny Point Classifier" << std::endl;

    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/stanford_bunny.stl");
    GeometricalEntityClassifier classifier(triangle_mesh);


    std::vector<PointType> rPoints{};
    rPoints.reserve(138600);
    for(double x = -24; x <= 85; x += 2){
        for(double y = -43; y <= 46; y += 2){
            for(double z = 5; z <= 115; z += 2){
                rPoints.push_back( {x, y, z} );
            }
        }
    }

    //auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<bool> result(rPoints.size(), false);
    int count = 0;
    for( auto& point : rPoints){
        if( classifier.IsInside(point) ){
            result[count] = true;
        }
        count++;
    }
    // auto end_time = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed_time = end_time - start_time;
    // std::cout << "Time new: " << elapsed_time.count() << std::endl;

    std::vector<bool> result_ref{};
    // Read reference results from file
    std::string line;
    std::ifstream myfile ("tibra/tests/cpp_tests/results/inside_outside_bunny.txt");
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
        result_ref.push_back( std::stoi(line) );
        }
        myfile.close();
    }

    // std::ofstream myfile;
    // myfile.open ("test.txt");
    // Compare results
    for( int i = 0; i < result.size(); ++i){
        BOOST_CHECK_EQUAL(result[i], result_ref[i]);
    }
    // myfile.close();
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing