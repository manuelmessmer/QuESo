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

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/io/io_utilities.h"
#include "queso/embedding/brep_operator.h"

#include "queso/tests/cpp_tests/global_config.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( ClipperTestSuite )

BOOST_AUTO_TEST_CASE(ClipCubeTest1) {
    QuESo_INFO << "Testing :: Test Clipper :: Clip Cylinder Test 1" << std::endl;

    //Read mesh from STL file
    TriangleMesh triangle_mesh{};

    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/cylinder.stl");

    Vector3d lower_bound = {0.0, 0.0, -0.1};
    Vector3d upper_bound = {2, 2, 1};

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    // Clip mesh.
    auto p_clipped_mesh = brep_operator.pClipTriangleMesh(lower_bound, upper_bound);

    double area = 0.0;
    for( IndexType i = 0; i < p_clipped_mesh->NumOfTriangles(); ++i){
        area += p_clipped_mesh->Area(i);
    }

    // Check surface area of clipped mesh.
    const double ref = 2.35619449;  // pi/2+pi/4 -> One quarte of the lateral surface + a quarter of the head face.
    QuESo_CHECK_NEAR( area, ref, 5e-4);
}

BOOST_AUTO_TEST_CASE(ClipCubeTest2) {
    QuESo_INFO << "Testing :: Test Clipper :: Clip Cylinder Test 2" << std::endl;

    //Read mesh from STL file
    TriangleMesh triangle_mesh{};

    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/cylinder.stl");

    Vector3d lower_bound = {0.0, 0.0, 0.0};
    Vector3d upper_bound = {2, 2, 1};

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    // Clip mesh.
    auto p_clipped_mesh = brep_operator.pClipTriangleMesh(lower_bound, upper_bound);

    double area = 0.0;
    for( IndexType i = 0; i < p_clipped_mesh->NumOfTriangles(); ++i){
        area += p_clipped_mesh->Area(i);
    }

    // Check surface area of clipped mesh.
    const double ref = 1.570796327;  // pi/2 -> One quarte of the lateral surface.
    QuESo_CHECK_NEAR( area, ref, 5e-4);
}

BOOST_AUTO_TEST_CASE(ClipCubeWithCavityTest) {
    QuESo_INFO << "Testing :: Test Clipper :: Clip Cube With Cavity Test" << std::endl;

    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/cube_with_cavity.stl");

    // Construct BRep_Operator
    BRepOperator brep_operator(triangle_mesh);

    const double delta_x = 0.15;
    const double delta_y = 0.15;
    const double delta_z = 0.15;

    std::vector<double> results_area{};
    for(double xx = -1.5001; xx <= 1.5; xx += delta_x){
        for(double yy = -1.5001; yy <= 1.5; yy += delta_y){
            for(double zz = -1.5001; zz <= 1.5; zz += delta_z){

                Vector3d lower_bound = {xx,yy,zz};
                Vector3d upper_bound = {xx+delta_x, yy+delta_y, zz+delta_z};
                // Clip mesh
                auto p_clipped_mesh = brep_operator.pClipTriangleMesh(lower_bound, upper_bound);
                // Compute area of clipped domain.
                double area = 0.0;
                for( IndexType i = 0; i < p_clipped_mesh->NumOfTriangles(); ++i){
                    area += p_clipped_mesh->Area(i);
                }
                results_area.push_back(area);
            }
        }
    }

    // Compare size
    QuESo_CHECK_EQUAL(results_area.size(), 9261 );

    // Compare reference areas of each clipped domain.
    std::ifstream file(base_dir + "/results/clipper_cube.txt");
    std::string line{};
    for( IndexType i = 0; i < results_area.size(); i++ ){
        std::getline(file, line);
        QuESo_CHECK_NEAR( std::stod(line), results_area[i], 1e-10);
    }
    file.close();
}

BOOST_AUTO_TEST_CASE(ClipElephantTest) {
    QuESo_INFO << "Testing :: Test Clipper :: Clip Elephant Test" << std::endl;

    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/elephant.stl");

    // Construct BRep_operator.
    BRepOperator brep_operator(triangle_mesh);

    const double delta_x = 0.05;
    const double delta_y = 0.05;
    const double delta_z = 0.05;

    std::vector<double> results_area{};
    for(double xx = -0.4; xx <= 0.4; xx += delta_x){
        for(double yy = -0.6; yy <= 0.6; yy += delta_y){
            for(double zz = -0.35; zz <= 0.35; zz += delta_z){

                Vector3d lower_bound = {xx,yy,zz};
                Vector3d upper_bound = {xx+delta_x, yy+delta_y, zz+delta_z};

                // Clip  Mesh.
                auto p_clipped_mesh = brep_operator.pClipTriangleMesh(lower_bound, upper_bound);

                // Compute area of clipped domain.
                double area = 0.0;
                for( IndexType i = 0; i < p_clipped_mesh->NumOfTriangles(); ++i){
                    area += p_clipped_mesh->Area(i);
                }
                results_area.push_back(area);
            }
        }
    }

    // Compare size
    QuESo_CHECK_EQUAL(results_area.size(), 6375 );

    // Compare reference areas of each clipped domain.
    std::ifstream file(base_dir + "/results/clipper_elephant.txt");
    std::string line{};
    for( IndexType i = 0; i < results_area.size(); i++ ){
        std::getline(file, line);
        QuESo_CHECK_NEAR( std::stod(line), results_area[i], 1e-10 );
    }
    file.close();
}


BOOST_AUTO_TEST_CASE(ClipBunnyTest) {
    QuESo_INFO << "Testing :: Test Clipper :: Clip Bunny Test" << std::endl;

    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/stanford_bunny.stl");

    // Construct BRep operator.
    BRepOperator brep_operator(triangle_mesh);
    const double delta_x = 5;
    const double delta_y = 5;
    const double delta_z = 5;

    std::vector<double> results_area{};
    for(double xx = -24; xx <= 85; xx += delta_x){
        for(double yy = -43; yy <= 46; yy += delta_y){
            for(double zz = 5; zz <= 115; zz += delta_z){
                Vector3d lower_bound = {xx,yy,zz};
                Vector3d upper_bound = {xx+delta_x, yy+delta_y, zz+delta_z};

                // Clip mesh.
                auto p_clipped_mesh = brep_operator.pClipTriangleMesh(lower_bound, upper_bound);

                // Compute area of clipped domain.
                double area = 0.0;
                for( IndexType i = 0; i < p_clipped_mesh->NumOfTriangles(); ++i){
                    area += p_clipped_mesh->Area(i);
                }
                results_area.push_back(area);
            }
        }
    }

    // Compare size
    QuESo_CHECK_EQUAL(results_area.size(), 9108 );

    // Compare reference areas of each clipped domain.
    std::ifstream file(base_dir + "/results/clipper_bunny.txt");
    std::string line{};
    for( IndexType i = 0; i < results_area.size(); i++ ){
        std::getline(file, line);
        QuESo_CHECK_NEAR( std::stod(line), results_area[i], 1e-10 );
    }
    file.close();
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso
