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
#include <doctest/doctest.h>
//// STL includes
#include <fstream>
#include <set>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/io/io_utilities.h"
#include "queso/quadrature/tensor_product.hpp"

#include "queso/tests/cpp/helpers/temporary_file.hpp"

namespace queso::testing {

    TEST_SUITE("IOTestSuite")
    {

        template<typename TMesh = TriangleMesh>
        TMesh CreateTestTriangleMesh()
        {
            TMesh mesh{};
            // Vertices
            mesh.AddVertex({ 0.0, 0.0, 0.0 });
            mesh.AddVertex({ 1.0, 0.0, 0.0 });
            mesh.AddVertex({ 0.0, 1.0, 0.0 });
            mesh.AddVertex({ 1.0, 1.0, 0.0 });
            // Triangles
            mesh.AddTriangle({ 0, 1, 2 }, { 0.0, 0.0, 1.0 });
            mesh.AddTriangle({ 1, 3, 2 }, { 0.0, 0.0, 1.0 });
            return mesh;
        }

        void TestTriangleMeshSTL(IO::EncodingType Encoding)
        {
            TemporaryFile tmp_file("temp_test_file.stl");

            auto mesh = CreateTestTriangleMesh();

            IO::WriteMeshToSTL(mesh, tmp_file.GetString(), Encoding);

            TriangleMesh mesh_read;
            IO::ReadMeshFromSTL(mesh_read, tmp_file.GetString());

            REQUIRE((mesh.NumOfTriangles()) == (mesh_read.NumOfTriangles()));
            REQUIRE((mesh.NumOfVertices()) == (mesh_read.NumOfVertices()));

            for (IndexType i = 0; i < mesh.NumOfTriangles(); ++i) {
                const auto tri_a = mesh.Triangle<WithoutNormals>(i);
                const auto tri_b = mesh_read.Triangle<WithoutNormals>(i);
                CheckPointNear(tri_a.P1, tri_b.P1, 1e-10);
                CheckPointNear(tri_a.P2, tri_b.P2, 1e-10);
                CheckPointNear(tri_a.P3, tri_b.P3, 1e-10);
            }
        }

        TEST_CASE("IOTestSuite - IOTriangle Mesh STLTest")
        {

            TestTriangleMeshSTL(IO::EncodingType::binary);
            TestTriangleMeshSTL(IO::EncodingType::ascii);
        }

        TEST_CASE("IOTestSuite - IODictionary JSONEscapes Strings Test")
        {
            TemporaryFile temporary_file("queso_dictionary.json");
            auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
            p_settings->SetValue(MainSettings::component_name, std::string("main"));
            p_settings->SetValue(MainSettings::input_filename, std::string("C:\\Users\\test\\\"quote\"\n\x01"));

            IO::WriteDictionaryToJSON(*p_settings, temporary_file.GetString());

            std::ifstream file(temporary_file.GetPath());
            const std::string contents{ std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>() };
            CHECK(contents.find("\"input_filename\"") != std::string::npos);
            CHECK(contents.find("C:\\\\Users\\\\test") != std::string::npos);
            CHECK(contents.find("\\\"quote\\\"") != std::string::npos);
            CHECK(contents.find("\\n\\u0001") != std::string::npos);
        }

        void CheckVTKBlocks(const std::filesystem::path& rFilename, const std::vector<std::string>& rBlockNames)
        {

            std::ifstream file(rFilename, std::ios::binary);
            REQUIRE(file);

            std::string file_contents((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

            std::set<std::string> blocks_found;
            for (const auto& block_name : rBlockNames) {
                if (file_contents.find(block_name) != std::string::npos) { blocks_found.insert(block_name); }
            }

            for (const auto& block_name : rBlockNames) { CHECK(blocks_found.find(block_name) != blocks_found.end()); }
        }

        void TestTriangleMeshVTK(IO::EncodingType Encoding)
        {
            TemporaryFile tmp_file("temp_test_file.stl");

            auto mesh = CreateTestTriangleMesh();

            IO::WriteMeshToVTK(mesh, tmp_file.GetString(), Encoding);
            CheckVTKBlocks(tmp_file.GetPath(), { "CELLS", "POINTS", "CELL_TYPES" });
        }

        TEST_CASE("IOTestSuite - IOTriangle Mesh VTKTest")
        {

            TestTriangleMeshVTK(IO::EncodingType::binary);
            TestTriangleMeshVTK(IO::EncodingType::ascii);
        }

        namespace {
            using IntegrationPointType = IntegrationPoint;
            using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
            using BackgroundGridType = BackgroundGrid<IntegrationPointType, BoundaryIntegrationPointType>;

            struct GaussBuilder
            {
                static constexpr BackgroundGridType::ElementFilter Builds =
                    BackgroundGridType::ElementFilter::untrimmed;
                Vector3i order;
                std::optional<BackgroundGridType::UntrimmedElementType>
                    Build(IndexType cell_index, const ElementBounds& bounds)
                {
                    BackgroundGridType::UntrimmedElementType element(cell_index + 1, bounds);
                    quadrature::tensor_product::Compute(
                        element, { .integration_order = order, .method = IntegrationMethod::gauss }
                    );
                    return element;
                }
            };
        }  // namespace

        Unique<BackgroundGridType> CreateTestBackgroundGrid()
        {
            auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
            auto& r_settings = *p_settings;

            auto& r_background_grid_settings = r_settings[MainSettings::background_grid_settings];
            r_background_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i({ 1, 1, 2 }));
            r_background_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
            r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ -24, -43, 5 });
            r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{ 85, 46, 115 });
            r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ -1.0, -1.0, -1.0 });
            r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 1.0, 1.0, 1.0 });
            r_background_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
            r_background_grid_settings.CheckRequired();

            Unique<BackgroundGridType> p_grid = MakeUnique<BackgroundGridType>(r_settings);

            GaussBuilder builder{ Vector3i({ 1, 2, 1 }) };

            IndexType number_elements = 2;
            for (IndexType i = 0; i < number_elements; ++i) { p_grid->MakeElement(builder, i); }
            p_grid->LockElements();

            return p_grid;
        }

        void TestElementsVTK(IO::EncodingType Encoding)
        {
            TemporaryFile tmp_file("temp_test_file.stl");

            auto p_grid = CreateTestBackgroundGrid();

            IO::WriteElementsToVTK(*p_grid, tmp_file.GetString(), Encoding);
            CheckVTKBlocks(tmp_file.GetPath(), { "CELLS", "POINTS", "CELL_TYPES" });
        }

        TEST_CASE("IOTestSuite - IOElement VTKTest")
        {

            TestElementsVTK(IO::EncodingType::binary);
            TestElementsVTK(IO::EncodingType::ascii);
        }

        void TestPointsVTK(IO::EncodingType Encoding)
        {
            TemporaryFile tmp_file("temp_test_file.stl");

            auto p_grid = CreateTestBackgroundGrid();

            IO::WritePointsToVTK(*p_grid, tmp_file.GetString(), Encoding);
            CheckVTKBlocks(tmp_file.GetPath(), { "CELLS", "POINTS", "CELL_TYPES", "POINT_DATA" });
        }

        TEST_CASE("IOTestSuite - IOPoints VTKTest")
        {

            TestPointsVTK(IO::EncodingType::binary);
            TestPointsVTK(IO::EncodingType::ascii);
        }

        using ConditionType = Condition<BackgroundGridType::ElementViewType>;

        Unique<ConditionType> CreateTestConditions()
        {
            auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
            auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");

            auto p_condition = MakeUnique<ConditionType>(*p_cond_settings, *p_cond_info);
            auto mesh = CreateTestTriangleMesh<TriangleMesh>();

            auto segment = ConditionType::ConditionSegmentType(1u, std::move(mesh));
            CHECK((segment.GetParentCellIndex()) == (IndexType{ 1 }));
            CHECK(!segment.IsInActiveElement());
            REQUIRE(!segment.GetParentElement().has_value());

            p_condition->AddSegment(std::move(segment));

            return p_condition;
        }

        void TestConditionSTL(IO::EncodingType Encoding)
        {
            TemporaryFile tmp_file("temp_test_file.stl");

            auto p_condition = CreateTestConditions();
            IO::WriteConditionToSTL(*p_condition, tmp_file.GetString(), Encoding);

            TriangleMesh mesh_read;
            IO::ReadMeshFromSTL(mesh_read, tmp_file.GetString());

            const auto& r_mesh = p_condition->GetSegments().begin()->GetTriangleMesh();

            REQUIRE((r_mesh.NumOfTriangles()) == (mesh_read.NumOfTriangles()));
            REQUIRE((r_mesh.NumOfVertices()) == (mesh_read.NumOfVertices()));

            for (IndexType i = 0; i < r_mesh.NumOfTriangles(); ++i) {
                const auto tri_a = r_mesh.Triangle<WithoutNormals>(i);
                const auto tri_b = mesh_read.Triangle<WithoutNormals>(i);
                CheckPointNear(tri_a.P1, tri_b.P1, 1e-10);
                CheckPointNear(tri_a.P2, tri_b.P2, 1e-10);
                CheckPointNear(tri_a.P3, tri_b.P3, 1e-10);
            }
        }

        TEST_CASE("IOTestSuite - IOCondition STLTest")
        {

            TestConditionSTL(IO::EncodingType::binary);
            TestConditionSTL(IO::EncodingType::ascii);
        }

    }  // TEST_SUITE

} // namespace queso::testing
