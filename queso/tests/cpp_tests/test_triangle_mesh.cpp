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
#include <cmath>
//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/clipped_triangle_mesh.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/brep_operator.h"
#include "queso/utilities/triangle_utilities.hpp"
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
    double surface_area = MeshUtilities::Area(triangle_mesh.View());
    double surface_area_omp = MeshUtilities::AreaOMP(triangle_mesh.View());

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
    double surface_area = MeshUtilities::Area(triangle_mesh.View());
    double surface_area_omp = MeshUtilities::AreaOMP(triangle_mesh.View());

    QuESo_CHECK_NEAR(surface_area, 69.11212872984862, 1e-10);
    QuESo_CHECK_NEAR(surface_area_omp, 69.11212872984862, 1e-10);
}

std::pair<double,Vector3d> ComputeAreaAndWeightedNormal(const TriangleMeshView& rTriangleMesh){
    double area = 0.0;
    Vector3d weighted_normal{0.0, 0.0, 0.0};

    rTriangleMesh.VisitEachTriangle<WithNormals>([&](const auto &triangle){
        area += TriangleUtilities::Area(triangle);
        const auto ips = TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPoint>(triangle, 1);
        for(const auto &ip : ips){
            const Vector3d n{triangle.Normal[0], triangle.Normal[1], triangle.Normal[2]};
            weighted_normal += (ip.Weight() * n);
        }
    });
    return std::make_pair(area, weighted_normal);
}

BOOST_AUTO_TEST_CASE(ClippedTriangleMeshBasicTest) {
    ClippedTriangleMesh mesh{};
    const IndexType v0 = mesh.NumOfVertices();
    mesh.AddVertex({0.0, 0.0, 0.0});
    const IndexType v1 = mesh.NumOfVertices();
    mesh.AddVertex({1.0, 0.0, 0.0});
    const IndexType v2 = mesh.NumOfVertices();
    mesh.AddVertex({0.0, 1.0, 0.0});
    mesh.AddTriangle({v0, v1, v2}, {0.0, 0.0, 1.0});
    mesh.AddEdgeOnPlane(SignedAxis::pos_z, v0, v1, 0);

    QuESo_CHECK_EQUAL(mesh.NumOfTriangles(), 1);
    QuESo_CHECK_EQUAL(mesh.NumberOfEdges(SignedAxis::pos_z), 1);

    IndexType edges = 0;
    for(const auto &edge : mesh.Edges(SignedAxis::pos_z)) {
        ++edges;
        QuESo_CHECK_NEAR(edge.mP1[0], 0.0, 1e-15);
        QuESo_CHECK_NEAR(edge.mP2[0], 1.0, 1e-15);
        QuESo_CHECK_NEAR(edge.mNormal[2], 1.0, 1e-15);
    }
    QuESo_CHECK_EQUAL(edges, 1);

    const auto view = mesh.MeshView();
    QuESo_CHECK_EQUAL(view.NumOfTriangles(), 1);
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
    auto init_values = ComputeAreaAndWeightedNormal(triangle_mesh.View());

    MeshUtilities::Refine(triangle_mesh, 20000);
    triangle_mesh.Check();
    const IndexType num_triangles = triangle_mesh.NumOfTriangles();
    QuESo_CHECK_GT(num_triangles, 19999);
    QuESo_CHECK_LT(num_triangles, 25001);
    auto new_values = ComputeAreaAndWeightedNormal(triangle_mesh.View());

    double weighted_normal_error = Math::Norm( (init_values.second - new_values.second) );
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
    auto init_values = ComputeAreaAndWeightedNormal(triangle_mesh.View());

    TriangleMesh new_mesh{};
    MeshUtilities::Append(new_mesh, triangle_mesh);
    new_mesh.Check();
    QuESo_CHECK_EQUAL(new_mesh.NumOfTriangles(), 112402);
    auto new_values = ComputeAreaAndWeightedNormal(new_mesh.View());

    double weighted_normal_error = Math::Norm( (init_values.second - new_values.second) );
    QuESo_CHECK_LT( weighted_normal_error, 1e-14);
    QuESo_CHECK_LT( std::abs(init_values.first - new_values.first)/init_values.first, 1e-10 );
}

BOOST_AUTO_TEST_CASE(TriangleMeshComputeVolumeBunnyTest) {
    QuESo_INFO << "Testing :: Test Triangle Mesh :: Test Compute Volume Bunny" << std::endl;
    TriangleMesh triangle_mesh{};
    // Read mesh from STL file
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/data/stanford_bunny.stl");
    double volume = MeshUtilities::Volume(triangle_mesh.View());
    double volume_omp = MeshUtilities::VolumeOMP(triangle_mesh.View());

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

    double volume = MeshUtilities::Volume(triangle_mesh.View());
    double volume_omp = MeshUtilities::VolumeOMP(triangle_mesh.View());

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

    double volume = MeshUtilities::Volume(triangle_mesh.View());
    double volume_omp = MeshUtilities::VolumeOMP(triangle_mesh.View());

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
                    const auto& r_mesh = p_trimmed_domain->GetTriangleMesh();
                    volume += MeshUtilities::Volume(r_mesh.View());
                }
                else if (status == IntersectionState::inside ){
                    const auto p_cube_mesh = MeshUtilities::pGetCuboid(lower_bound, upper_bound);
                    volume += MeshUtilities::Volume(p_cube_mesh->View());
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

    double volume = MeshUtilities::Volume(triangle_mesh.View());
    double volume_omp = MeshUtilities::VolumeOMP(triangle_mesh.View());

    const double volume_ref = 22.81560787501277;
    QuESo_CHECK_LT( std::abs(volume - volume_ref) / volume_ref, 1e-9);
    QuESo_CHECK_LT( std::abs(volume_omp - volume_ref) / volume_ref, 1e-9);
}

void CheckTriangleMeshViewTraversalBranches(const TriangleMeshView& rView)
{
    QuESo_CHECK_EQUAL(rView.NumOfTriangles(), 4);

    double area_visit_all = 0.0;
    rView.Visit<WithNormals>([&](const auto& triangles) {
        for (const auto& triangle : triangles) {
            area_visit_all += TriangleUtilities::Area(triangle);
        }
    });

    IndexType num_without_normals = 0;
    double area_each_without_normals = 0.0;
    rView.VisitEachTriangle<WithoutNormals>([&](const auto& triangle) {
        ++num_without_normals;
        area_each_without_normals += TriangleUtilities::Area(triangle);
    });

    const auto subset_indices = std::views::iota(IndexType{0}, IndexType{2});
    double area_visit_subset = 0.0;
    rView.Visit<WithNormals>(subset_indices, [&](const auto& triangles) {
        for (const auto& triangle : triangles) {
            area_visit_subset += TriangleUtilities::Area(triangle);
        }
    });

    IndexType visited_until_stop = 0;
    rView.VisitEachTriangle<WithNormals>([&](const auto&) {
        ++visited_until_stop;
        return TriangleMeshView::VisitToken::stop_loop;
    });

    QuESo_CHECK_NEAR(area_visit_all, 2.3660254037844384, 1e-14);
    QuESo_CHECK_NEAR(area_each_without_normals, 2.3660254037844384, 1e-14);
    QuESo_CHECK_NEAR(area_visit_subset, 1.0, 1e-14);
    QuESo_CHECK_EQUAL(num_without_normals, 4);
    QuESo_CHECK_EQUAL(visited_until_stop, 1);
}

BOOST_AUTO_TEST_CASE(TriangleMeshViewNativeTraversalTest) {
    TriangleMesh native_mesh{};
    const IndexType n0 = native_mesh.AddVertex({0.0, 0.0, 0.0});
    const IndexType n1 = native_mesh.AddVertex({1.0, 0.0, 0.0});
    const IndexType n2 = native_mesh.AddVertex({0.0, 1.0, 0.0});
    const IndexType n3 = native_mesh.AddVertex({0.0, 0.0, 1.0});
    native_mesh.AddTriangle({n0, n1, n2}, {0.0, 0.0, 1.0});
    native_mesh.AddTriangle({n0, n1, n3}, {0.0, -1.0, 0.0});
    native_mesh.AddTriangle({n0, n2, n3}, {1.0, 0.0, 0.0});
    native_mesh.AddTriangle({n1, n2, n3}, Vector3d{1.0, 1.0, 1.0} / std::sqrt(3.0));

    CheckTriangleMeshViewTraversalBranches(native_mesh.View());
}

BOOST_AUTO_TEST_CASE(TriangleMeshViewForeignTraversalTest) {
    struct ForeignMesh {
        std::vector<Vector3d> vertices;
        std::vector<Vector3i> triangles;
    };

    ForeignMesh mesh{
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}},
        {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}
    };

    TriangleMeshView view(
        mesh,
        [](const ForeignMesh &m){ return m.triangles.size(); },
        [](const ForeignMesh &m, IndexType i){
            const auto tri = m.triangles[i];
            return TriangleProxy<WithoutNormals>{m.vertices[tri[0]], m.vertices[tri[1]], m.vertices[tri[2]]};
        }
    );

    CheckTriangleMeshViewTraversalBranches(view);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

