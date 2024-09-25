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
//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/containers/background_grid.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/brep_operator.h"
#include "queso/io/io_utilities.h"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( BoundingBoxOfTrimmedDomainTestSuite )

BOOST_AUTO_TEST_CASE(CylinderBoundingBoxOfTrimmedDomainTest) {

    QuESo_INFO << "Testing :: Test Bounding Box Of Trimmed Domain:: Cylinder" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cylinder.stl");

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double delta_x = 0.50;
    const double delta_y = 0.50;
    const double delta_z = 0.50;

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    std::vector<double> results;
    results.reserve(344);
    for(double x = -1.501; x <= 1.5; x += delta_x){
        for(double y = -1.501; y <= 1.5; y += delta_y){
            for(double z = -1.01; z <= 12; z += delta_z){
                Vector3d lower_bound = {x, y, z};
                Vector3d upper_bound = {x+delta_x, y+delta_y, z+delta_z};
                auto status = brep_operator.GetIntersectionState(lower_bound, upper_bound);
                if( status == IntersectionStatus::Trimmed){
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound, min_vol_ratio, min_num_triangles);
                    auto bounding_box = p_trimmed_domain->GetBoundingBoxOfTrimmedDomain();
                    const double norm = Math::Norm( Math::Subtract(bounding_box.second, bounding_box.first) );
                    results.push_back( norm );
                }
            }
        }
    }

    QuESo_CHECK_EQUAL(results.size(), 344);
    std::ifstream myfile("queso/tests/cpp_tests/results/aabb_trimmed_domain_cylinder.txt");
    std::string line;
    for( IndexType i = 0; i < results.size(); ++i){
        getline (myfile, line);
        double ref_value = std::stod(line);
        double error = std::abs(results[i] - ref_value)/ref_value;
        QuESo_CHECK_LT( error, 1e-10 );
    }
    myfile.close();
}

BOOST_AUTO_TEST_CASE(CubeBoundingBoxOfTrimmedDomainTest) {
    QuESo_INFO << "Testing :: Test Bounding Box Of Trimmed Domain:: Cube" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cube_with_cavity.stl");

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double delta_x = 0.3;
    const double delta_y = 0.3;
    const double delta_z = 0.3;

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    std::vector<double> results{};
    results.reserve(826);
    for(double x = -1.501; x <= 1.5; x += delta_x){
        for(double y = -1.501; y <= 1.5; y += delta_y){
            for(double z = -1.501; z <= 1.5; z += delta_z){
                Vector3d lower_bound = {x, y, z};
                Vector3d upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto status = brep_operator.GetIntersectionState(lower_bound, upper_bound);
                if( status == IntersectionStatus::Trimmed){
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound, min_vol_ratio, min_num_triangles);
                    auto bounding_box = p_trimmed_domain->GetBoundingBoxOfTrimmedDomain();
                    const double norm = Math::Norm( Math::Subtract(bounding_box.second, bounding_box.first) );
                    results.push_back( norm );
                }
            }
        }
    }

    QuESo_CHECK_EQUAL(results.size(), 826);
    std::ifstream myfile("queso/tests/cpp_tests/results/aabb_trimmed_domain_cube.txt");
    std::string line;
    for( IndexType i = 0; i < results.size(); ++i){
        getline (myfile, line);
        double ref_value = std::stod(line);
        double error = std::abs(results[i] - ref_value)/ref_value;
        QuESo_CHECK_LT( error, 1e-10 );
    }
    myfile.close();

}

BOOST_AUTO_TEST_CASE(ElephantBoundingBoxOfTrimmedDomainTest) {
    QuESo_INFO << "Testing :: Test Bounding Box Of Trimmed Domain:: Elephant" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/elephant.stl");

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double delta_x = 0.1;
    const double delta_y = 0.1;
    const double delta_z = 0.1;

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    std::vector<double> results{};
    results.reserve(166);
    for(double x = -0.4; x <= 0.4; x += delta_x){
        for(double y = -0.6; y <= 0.6; y += delta_y){
            for(double z = -0.35; z <= 0.35; z += delta_x){
                Vector3d lower_bound = {x, y, z};
                Vector3d upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto status = brep_operator.GetIntersectionState(lower_bound, upper_bound);
                if( status == IntersectionStatus::Trimmed){
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound, min_vol_ratio, min_num_triangles);
                    auto bounding_box = p_trimmed_domain->GetBoundingBoxOfTrimmedDomain();
                    const double norm = Math::Norm( Math::Subtract(bounding_box.second, bounding_box.first) );
                    results.push_back( norm );
                }
            }
        }
    }

    QuESo_CHECK_EQUAL(results.size(), 166);
    std::ifstream myfile("queso/tests/cpp_tests/results/aabb_trimmed_domain_elephant.txt");
    std::string line;
    for( IndexType i = 0; i < results.size(); ++i){
        getline (myfile, line);
        double ref_value = std::stod(line);
        double error = std::abs(results[i] - ref_value)/ref_value;
        QuESo_CHECK_LT( error, 1e-10 );
    }
    myfile.close();

}

BOOST_AUTO_TEST_CASE(BunnyBoundingBoxOfTrimmedDomainTest) {
    QuESo_INFO << "Testing :: Test Bounding Box Of Trimmed Domain:: Bunny" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/stanford_bunny.stl");

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double delta_x = 10;
    const double delta_y = 10;
    const double delta_z = 10;

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    std::vector<double> results{};
    results.reserve(381);
    for(double x = -24; x <= 85; x += delta_x){
        for(double y = -43; y <= 46; y += delta_y){
            for(double z = 5; z <= 115; z += delta_z){
                Vector3d lower_bound = {x, y, z};
                Vector3d upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto status = brep_operator.GetIntersectionState(lower_bound, upper_bound);
                if( status == IntersectionStatus::Trimmed){
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound, min_vol_ratio, min_num_triangles);
                    auto bounding_box = p_trimmed_domain->GetBoundingBoxOfTrimmedDomain();
                    const double norm = Math::Norm( Math::Subtract(bounding_box.second, bounding_box.first) );
                    results.push_back( norm );
                }
            }
        }
    }
    QuESo_CHECK_EQUAL(results.size(), 381);
    std::ifstream myfile("queso/tests/cpp_tests/results/aabb_trimmed_domain_bunny.txt");
    std::string line;
    for( IndexType i = 0; i < results.size(); ++i){
        getline (myfile, line);
        double ref_value = std::stod(line);
        double error = std::abs(results[i] - ref_value)/ref_value;
        QuESo_CHECK_LT( error, 1e-10 );
    }
    myfile.close();
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso