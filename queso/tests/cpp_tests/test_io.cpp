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
//// STL includes
#include <set>

//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/io/io_utilities.h"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/quadrature/single_element.hpp"

#include "queso/tests/cpp_tests/helper/temporary_file.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( IOTestSuite )

template<typename TMesh = TriangleMesh>
TMesh CreateTestTriangleMesh() {
    TMesh mesh{};
    // Vertices
    mesh.AddVertex({0.0, 0.0, 0.0});
    mesh.AddVertex({1.0, 0.0, 0.0});
    mesh.AddVertex({0.0, 1.0, 0.0});
    mesh.AddVertex({1.0, 1.0, 0.0});
    // Triangles
    mesh.AddTriangle({0, 1, 2}, {0.0, 0.0, 1.0});
    mesh.AddTriangle({1, 3, 2}, {0.0, 0.0, 1.0});
    return mesh;
}

void TestTriangleMeshSTL(IO::EncodingType Encoding) {
    TemporaryFile tmp_file("temp_test_file.stl");

    auto mesh = CreateTestTriangleMesh();

    IO::WriteMeshToSTL(mesh, tmp_file.GetString(), Encoding);

    TriangleMesh mesh_read;
    IO::ReadMeshFromSTL(mesh_read, tmp_file.GetString());

    QuESo_CHECK_EQUAL(mesh.NumOfTriangles(), mesh_read.NumOfTriangles());
    QuESo_CHECK_EQUAL(mesh.NumOfVertices(), mesh_read.NumOfVertices());

    for(IndexType i = 0; i < mesh.NumOfTriangles(); ++i) {
        const auto tri_a = mesh.Triangle<WithoutNormals>(i);
        const auto tri_b = mesh_read.Triangle<WithoutNormals>(i);
        QuESo_CHECK_POINT_NEAR( tri_a.P1, tri_b.P1, 1e-10 );
        QuESo_CHECK_POINT_NEAR( tri_a.P2, tri_b.P2, 1e-10 );
        QuESo_CHECK_POINT_NEAR( tri_a.P3, tri_b.P3, 1e-10 );
    }
}

BOOST_AUTO_TEST_CASE(IOTriangleMeshSTLTest) {
    QuESo_INFO << "Testing :: Test IO Triangle Mesh :: Read / Write STL" << std::endl;
    TestTriangleMeshSTL(IO::EncodingType::binary);
    TestTriangleMeshSTL(IO::EncodingType::ascii);
}

void CheckVTKBlocks(const std::filesystem::path& rFilename,
                    const std::vector<std::string>& rBlockNames) {

    std::ifstream file(rFilename, std::ios::binary);
    BOOST_REQUIRE(file);

    std::string file_contents((std::istreambuf_iterator<char>(file)),
                               std::istreambuf_iterator<char>());

    std::set<std::string> blocks_found;
    for (const auto& block_name : rBlockNames) {
        if (file_contents.find(block_name) != std::string::npos) {
            blocks_found.insert(block_name);
        }
    }

    for (const auto& block_name : rBlockNames) {
        QuESo_CHECK(blocks_found.find(block_name) != blocks_found.end());
    }
}

void TestTriangleMeshVTK(IO::EncodingType Encoding) {
    TemporaryFile tmp_file("temp_test_file.stl");

    auto mesh = CreateTestTriangleMesh();

    IO::WriteMeshToVTK(mesh, tmp_file.GetString(), Encoding);
    CheckVTKBlocks(tmp_file.GetPath(), {"CELLS", "POINTS", "CELL_TYPES"});
}

BOOST_AUTO_TEST_CASE(IOTriangleMeshVTKTest) {
    QuESo_INFO << "Testing :: Test IO Triangle Mesh :: Write VTK" << std::endl;
    TestTriangleMeshVTK(IO::EncodingType::binary);
    TestTriangleMeshVTK(IO::EncodingType::ascii);
}

namespace {
	using IntegrationPointType = IntegrationPoint;
	using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
	using BackgroundGridType = BackgroundGrid<IntegrationPointType, BoundaryIntegrationPointType>;

	struct GaussBuilder { 
		static constexpr BackgroundGridType::ElementFilter Builds = BackgroundGridType::ElementFilter::untrimmed;
		Vector3i order;
		std::optional<BackgroundGridType::UntrimmedElementType> Build(IndexType id, const ElementBounds& bounds) {
			BackgroundGridType::UntrimmedElementType element(id, bounds);
			QuadratureSingleElement<BackgroundGridType::UntrimmedElementType>::AssembleIPs(
				element, order, IntegrationMethod::gauss
			);
			return element;
		}
	};
}

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
    r_background_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{2, 2, 2});
	r_background_grid_settings.CheckRequired();

    Unique<BackgroundGridType> p_grid = MakeUnique<BackgroundGridType>(r_settings);

    GaussBuilder builder{ Vector3i({1, 2, 1}) };

    IndexType number_elements = 2;
    for( IndexType i = 1; i <= number_elements; ++i){
        const auto bounds_xyz = MakeBox({0.0, 0.0, 0.0}, {0.1, 0.1, 0.1});
        const auto bounds_uvw = MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
        p_grid->MakeElement(builder, i, ElementBounds{bounds_xyz, bounds_uvw});
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
    CheckVTKBlocks(tmp_file.GetPath(), {"CELLS", "POINTS", "CELL_TYPES", "POINT_DATA"});
}

BOOST_AUTO_TEST_CASE(IOPointsVTKTest) {
    QuESo_INFO << "Testing :: Test IO Points :: Write VTK" << std::endl;
    TestPointsVTK(IO::EncodingType::binary);
    TestPointsVTK(IO::EncodingType::ascii);
}

using ConditionType = Condition<BackgroundGridType::ElementViewType>;

Unique<ConditionType> CreateTestConditions(){
    auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
    auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");

    auto p_condition = MakeUnique<ConditionType>(*p_cond_settings, *p_cond_info);
    auto mesh = CreateTestTriangleMesh<ClippedTriangleMesh>();

    auto segment = ConditionType::ConditionSegmentType(1u, std::move(mesh));

    p_condition->AddSegment(std::move(segment));

    return p_condition;
}

void TestConditionSTL(IO::EncodingType Encoding) {
    TemporaryFile tmp_file("temp_test_file.stl");

    auto p_condition = CreateTestConditions();
    IO::WriteConditionToSTL(*p_condition, tmp_file.GetString(), Encoding);

    TriangleMesh mesh_read;
    IO::ReadMeshFromSTL(mesh_read, tmp_file.GetString());

    const auto& r_mesh = p_condition->GetSegments().begin()->GetTriangleMesh();

    QuESo_CHECK_EQUAL(r_mesh.NumOfTriangles(), mesh_read.NumOfTriangles());
    QuESo_CHECK_EQUAL(r_mesh.NumOfVertices(), mesh_read.NumOfVertices());

    for(IndexType i = 0; i < r_mesh.NumOfTriangles(); ++i) {
        const auto tri_a = r_mesh.Triangle<WithoutNormals>(i);
        const auto tri_b = mesh_read.Triangle<WithoutNormals>(i);
        QuESo_CHECK_POINT_NEAR( tri_a.P1, tri_b.P1, 1e-10 );
        QuESo_CHECK_POINT_NEAR( tri_a.P2, tri_b.P2, 1e-10 );
        QuESo_CHECK_POINT_NEAR( tri_a.P3, tri_b.P3, 1e-10 );
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
