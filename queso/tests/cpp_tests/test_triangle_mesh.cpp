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
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/brep_operator.h"
#include "queso/utilities/mesh_utilities.h"
#include "queso/io/io_utilities.h"

#include "queso/tests/cpp_tests/global_config.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( TriangleMeshTestSuite )

BOOST_AUTO_TEST_CASE(TriangleMeshIOBindaryTest) {
    QuESo_INFO << "Testing :: Test Triangle Mesh :: Test IO Binary" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/cylinder.stl");

    // Make basic check
    triangle_mesh.Check();

    QuESo_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 888);
    // Check surface area
    double surface_area = MeshUtilities::Area(triangle_mesh);
    double surface_area_omp = MeshUtilities::AreaOMP(triangle_mesh);

    QuESo_CHECK_NEAR(surface_area, 69.11212872984862, 1e-10);
    QuESo_CHECK_NEAR(surface_area_omp, 69.11212872984862, 1e-10);
}

BOOST_AUTO_TEST_CASE(TriangleMeshIOAsciiTest) {
    QuESo_INFO << "Testing :: Test Triangle Mesh :: Test IO Ascii" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/cylinder_ascii.stl");

    // Make basic check
    triangle_mesh.Check();

    QuESo_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 888);
    // Check surface area
    double surface_area = MeshUtilities::Area(triangle_mesh);
    double surface_area_omp = MeshUtilities::AreaOMP(triangle_mesh);

    QuESo_CHECK_NEAR(surface_area, 69.11212872984862, 1e-10);
    QuESo_CHECK_NEAR(surface_area_omp, 69.11212872984862, 1e-10);
}

std::pair<double,Vector3d> ComputeAreaAndWeightedNormal(const TriangleMesh& rTriangleMesh){
    double area = 0.0;
    Vector3d weighted_normal{0.0, 0.0, 0.0};

    for( IndexType triangle_id = 0; triangle_id < rTriangleMesh.NumOfTriangles(); ++triangle_id){
        area += rTriangleMesh.Area(triangle_id);
        auto p_ips = rTriangleMesh.pGetIPsGlobal<BoundaryIntegrationPoint>(triangle_id, 1);
        for( auto& ip : (*p_ips) ){
            Math::AddSelf(weighted_normal, Math::Mult(ip.Weight(), rTriangleMesh.Normal(triangle_id)) );
        }
    }
    return std::make_pair(area, weighted_normal);
}

BOOST_AUTO_TEST_CASE(TriangleMeshRefineTest) {
    QuESo_INFO << "Testing :: Test Triangle Mesh :: Test Refine" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/elephant.stl");

    // Make basic check
    triangle_mesh.Check();
    QuESo_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 5558);
    auto init_values = ComputeAreaAndWeightedNormal(triangle_mesh);

    MeshUtilities::Refine(triangle_mesh, 20000);
    triangle_mesh.Check();
    QuESo_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 27911);
    auto new_values = ComputeAreaAndWeightedNormal(triangle_mesh);

    double weighted_normal_error = Math::Norm( Math::Subtract(init_values.second, new_values.second) );
    QuESo_CHECK_LT( weighted_normal_error, 1e-14);
    QuESo_CHECK_LT( std::abs(init_values.first - new_values.first)/init_values.first, 1e-10 );
}

BOOST_AUTO_TEST_CASE(TriangleMeshAppendTest) {
    QuESo_INFO << "Testing :: Test Triangle Mesh :: Test Append" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/stanford_bunny.stl");

    // Make basic check
    triangle_mesh.Check();
    QuESo_CHECK_EQUAL(triangle_mesh.NumOfTriangles(), 112402);
    auto init_values = ComputeAreaAndWeightedNormal(triangle_mesh);

    TriangleMesh new_mesh{};
    MeshUtilities::Append(new_mesh, triangle_mesh);
    new_mesh.Check();
    QuESo_CHECK_EQUAL(new_mesh.NumOfTriangles(), 112402);
    auto new_values = ComputeAreaAndWeightedNormal(new_mesh);

    double weighted_normal_error = Math::Norm( Math::Subtract(init_values.second, new_values.second) );
    QuESo_CHECK_LT( weighted_normal_error, 1e-14);
    QuESo_CHECK_LT( std::abs(init_values.first - new_values.first)/init_values.first, 1e-10 );
}

BOOST_AUTO_TEST_CASE(TriangleMeshComputeVolumeBunnyTest) {
    QuESo_INFO << "Testing :: Test Triangle Mesh :: Test Compute Volume Bunny" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/stanford_bunny.stl");
    double volume = MeshUtilities::Volume(triangle_mesh);
    double volume_omp = MeshUtilities::VolumeOMP(triangle_mesh);

    const double volume_ref = 279628.2991519215;
    QuESo_CHECK_LT( std::abs(volume - volume_ref) / volume_ref, 1e-9);
    QuESo_CHECK_LT( std::abs(volume_omp - volume_ref) / volume_ref, 1e-9);
}

BOOST_AUTO_TEST_CASE(TriangleMeshComputeCylinderTest) {
    QuESo_INFO << "Testing :: Test Triangle Mesh :: Test Compute Volume Cylinder" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/cylinder.stl");

    double volume = MeshUtilities::Volume(triangle_mesh);
    double volume_omp = MeshUtilities::VolumeOMP(triangle_mesh);

    const double volume_ref = 31.41176999044123;
    QuESo_CHECK_LT( std::abs(volume - volume_ref) / volume_ref, 1e-9);
    QuESo_CHECK_LT( std::abs(volume_omp - volume_ref) / volume_ref, 1e-9);
}

BOOST_AUTO_TEST_CASE(TriangleMeshComputeElephantTest) {
    QuESo_INFO << "Testing :: Test Triangle Mesh :: Test Compute Volume Elephant" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/elephant.stl");

    double volume = MeshUtilities::Volume(triangle_mesh);
    double volume_omp = MeshUtilities::VolumeOMP(triangle_mesh);

    const double volume_ref = 0.04620123478735502;
    QuESo_CHECK_LT( std::abs(volume - volume_ref) / volume_ref, 1e-8);
    QuESo_CHECK_LT( std::abs(volume_omp - volume_ref) / volume_ref, 1e-8);
}

BOOST_AUTO_TEST_CASE(TriangleMeshComputeElephant2Test) {
    QuESo_INFO << "Testing :: Test Triangle Mesh :: Test Compute Volume Elephant Splitted" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/elephant.stl");

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    const double delta_x = 0.1;
    const double delta_y = 0.1;
    const double delta_z = 0.1;

    std::vector<double> results{};
    results.reserve(166);
    double volume = 0.0;
    // Compute weight of each individual element.
    for(double x = -0.4; x <= 0.4; x += delta_x){
        for(double y = -0.6; y <= 0.6; y += delta_y){
            for(double z = -0.35; z <= 0.35; z += delta_x){
                Vector3d lower_bound = {x, y, z};
                Vector3d upper_bound = {x+delta_x, y+delta_y, z+delta_z};
                auto status = brep_operator.GetIntersectionState(lower_bound, upper_bound);
                if( status == IntersectionState::trimmed){
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound, upper_bound, min_vol_ratio, min_num_triangles);
                    const auto p_boundary_ips = p_trimmed_domain->pGetBoundaryIps<BoundaryIntegrationPoint>();
                    const auto& r_mesh = p_trimmed_domain->GetTriangleMesh();
                    volume += MeshUtilities::Volume(r_mesh);
                }
                else if (status == IntersectionState::inside ){
                    const auto p_cube_mesh = MeshUtilities::pGetCuboid(lower_bound, upper_bound);
                    volume += MeshUtilities::Volume(*p_cube_mesh);
                }
            }
        }
    }
    const double volume_ref = 0.04620123478735502;
    QuESo_CHECK_LT( std::abs(volume - volume_ref) / volume_ref, 1e-8);
}

BOOST_AUTO_TEST_CASE(TriangleMeshComputeCubeTest) {
    QuESo_INFO << "Testing :: Test Triangle Mesh :: Test Compute Volume Cube" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/cube_with_cavity.stl");

    double volume = MeshUtilities::Volume(triangle_mesh);
    double volume_omp = MeshUtilities::VolumeOMP(triangle_mesh);

    const double volume_ref = 22.81560787501277;
    QuESo_CHECK_LT( std::abs(volume - volume_ref) / volume_ref, 1e-9);
    QuESo_CHECK_LT( std::abs(volume_omp - volume_ref) / volume_ref, 1e-9);
}
BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso