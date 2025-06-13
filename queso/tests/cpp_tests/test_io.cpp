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
#include "queso/includes/dictionary_factory.hpp"
#include "queso/io/io_utilities.h"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/quadrature/single_element.hpp"

#include "queso/tests/cpp_tests/helper/temporary_file.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( IOTestSuite )

Unique<TriangleMeshInterface> CreateTestTriangleMesh() {
    auto p_mesh = MakeUnique<TriangleMesh>();
    // Vertices
    p_mesh->AddVertex({0.0, 0.0, 0.0});
    p_mesh->AddVertex({1.0, 0.0, 0.0});
    p_mesh->AddVertex({0.0, 1.0, 0.0});
    p_mesh->AddVertex({1.0, 1.0, 0.0});
    // Triangles
    p_mesh->AddTriangle({0, 1, 2});
    p_mesh->AddTriangle({1, 3, 2});
    // Normals
    p_mesh->AddNormal({0.0, 0.0, 1.0});
    p_mesh->AddNormal({0.0, 0.0, 1.0});
    return p_mesh;
}

void TestTriangleMeshSTL(IO::EncodingType Encoding) {
    TemporaryFile tmp_file("temp_test_file.stl");

    auto p_mesh = CreateTestTriangleMesh();

    IO::WriteMeshToSTL(*p_mesh, tmp_file.GetString(), Encoding);

    TriangleMesh mesh_read;
    IO::ReadMeshFromSTL(mesh_read, tmp_file.GetString());

    QuESo_CHECK_EQUAL(p_mesh->NumOfTriangles(), mesh_read.NumOfTriangles());
    QuESo_CHECK_EQUAL(p_mesh->NumOfVertices(), mesh_read.NumOfVertices());

    for(IndexType i = 0; i < p_mesh->NumOfTriangles(); ++i) {
        QuESo_CHECK_POINT_NEAR( p_mesh->P1(i), mesh_read.P1(i), 1e-10 );
        QuESo_CHECK_POINT_NEAR( p_mesh->P2(i), mesh_read.P2(i), 1e-10 );
        QuESo_CHECK_POINT_NEAR( p_mesh->P3(i), mesh_read.P3(i), 1e-10 );
    }
}

BOOST_AUTO_TEST_CASE(IOTriangleMeshSTLTest) {
    QuESo_INFO << "Testing :: Test IO Triangle Mesh :: Read / Write STL" << std::endl;
    TestTriangleMeshSTL(IO::EncodingType::binary);
    TestTriangleMeshSTL(IO::EncodingType::ascii);
}

void CheckVTKBlocks(const std::filesystem::path& rFilename, const std::vector<std::string>& rBlockNames) {
    std::ifstream file(rFilename);
    BOOST_REQUIRE(file.is_open());

    std::set<std::string> blocks_found;
    std::string line;
    while (std::getline(file, line)) {
        for(const auto& block_name : rBlockNames) {
            if( line.find(block_name) != std::string::npos ) {
                blocks_found.insert(block_name);
            }
        }
    }
    for(const auto& block_name : rBlockNames) {
        QuESo_CHECK( blocks_found.find(block_name) != blocks_found.end() );
    }
}

void TestTriangleMeshVTK(IO::EncodingType Encoding) {
    TemporaryFile tmp_file("temp_test_file.stl");

    auto p_mesh = CreateTestTriangleMesh();

    IO::WriteMeshToVTK(*p_mesh, tmp_file.GetString(), Encoding);
    CheckVTKBlocks(tmp_file.GetPath(), {"CELLS", "POINTS", "CELL_TYPES"});
}

BOOST_AUTO_TEST_CASE(IOTriangleMeshVTKTest) {
    QuESo_INFO << "Testing :: Test IO Triangle Mesh :: Write VTK" << std::endl;
    TestTriangleMeshSTL(IO::EncodingType::binary);
    TestTriangleMeshSTL(IO::EncodingType::ascii);
}

using IntegrationPointType = IntegrationPoint;
using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
using ElementType = Element<IntegrationPointType, BoundaryIntegrationPointType>;
using BackgroundGridType = BackgroundGrid<ElementType>;

Unique<BackgroundGridType> CreateTestBackgroundGrid(){
    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_background_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_background_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i({1,1,2}));
    r_background_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-24, -43, 5});
    r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{85, 46, 115});
    r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{-1.0, -1-0, 1.0});
    r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 1.0, 1.0});

    Unique<BackgroundGridType> p_grid = MakeUnique<BackgroundGridType>(r_settings);

    IndexType number_elements = 2;
    for( IndexType i = 1; i <= number_elements; ++i){
        PointType tmp_point_A = {0.0, 0.0, 0.0};
        PointType tmp_point_B = {0.1, 0.1, 0.1};
        Unique<ElementType> p_new_element = MakeUnique<ElementType>(i, MakeBox(tmp_point_A, tmp_point_B),
                                                                       MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}) );

        QuadratureSingleElement<ElementType>::AssembleIPs(*p_new_element, Vector3i({1, 2, 1}), IntegrationMethod::gauss);
        p_grid->AddElement(std::move(p_new_element));
    }

    return p_grid;
}

void TestElementsVTK(IO::EncodingType Encoding) {
    TemporaryFile tmp_file("temp_test_file.stl");

    auto p_grid = CreateTestBackgroundGrid();

    IO::WriteElementsToVTK(*p_grid, tmp_file.GetString(), Encoding);
    CheckVTKBlocks(tmp_file.GetPath(), {"CELLS", "POINTS", "CELL_TYPES"});
}

BOOST_AUTO_TEST_CASE(IOElementVTKTest) {
    QuESo_INFO << "Testing :: Test IO Elements :: Write VTK" << std::endl;
    TestElementsVTK(IO::EncodingType::binary);
    TestElementsVTK(IO::EncodingType::ascii);
}

void TestPointsVTK(IO::EncodingType Encoding) {
    TemporaryFile tmp_file("temp_test_file.stl");

    auto p_grid = CreateTestBackgroundGrid();

    IO::WritePointsToVTK(*p_grid, tmp_file.GetString(), Encoding);
    CheckVTKBlocks(tmp_file.GetPath(), {"CELLS", "POINTS", "CELL_TYPES"});
}

BOOST_AUTO_TEST_CASE(IOPointsVTKTest) {
    QuESo_INFO << "Testing :: Test IO Points :: Write VTK" << std::endl;
    TestPointsVTK(IO::EncodingType::binary);
    TestPointsVTK(IO::EncodingType::ascii);
}

using ConditionType = Condition<ElementType>;

Unique<ConditionType> CreateTestConditions(){
    auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
    auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");

    auto p_condition = MakeUnique<ConditionType>(*p_cond_settings, *p_cond_info);
    auto p_mesh = CreateTestTriangleMesh();

    auto p_segment = MakeUnique<ConditionType::ConditionSegmentType>(1u, p_mesh);

    p_condition->AddSegment(p_segment);

    return p_condition;
}

void TestConditionSTL(IO::EncodingType Encoding) {
    TemporaryFile tmp_file("temp_test_file.stl");

    auto p_condition = CreateTestConditions();
    IO::WriteConditionToSTL(*p_condition, tmp_file.GetString(), Encoding);

    TriangleMesh mesh_read;
    IO::ReadMeshFromSTL(mesh_read, tmp_file.GetString());

    const auto& r_mesh = p_condition->SegmentsBegin()->GetTriangleMesh();

    QuESo_CHECK_EQUAL(r_mesh.NumOfTriangles(), mesh_read.NumOfTriangles());
    QuESo_CHECK_EQUAL(r_mesh.NumOfVertices(), mesh_read.NumOfVertices());

    for(IndexType i = 0; i < r_mesh.NumOfTriangles(); ++i) {
        QuESo_CHECK_POINT_NEAR( r_mesh.P1(i), mesh_read.P1(i), 1e-10 );
        QuESo_CHECK_POINT_NEAR( r_mesh.P2(i), mesh_read.P2(i), 1e-10 );
        QuESo_CHECK_POINT_NEAR( r_mesh.P3(i), mesh_read.P3(i), 1e-10 );
    }
}

BOOST_AUTO_TEST_CASE(IOConditionSTLTest) {
    QuESo_INFO << "Testing :: Test IO Condition :: Write STL" << std::endl;
    TestConditionSTL(IO::EncodingType::binary);
    TestConditionSTL(IO::EncodingType::ascii);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso