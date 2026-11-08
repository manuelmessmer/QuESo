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
//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/domain_mesh_embedder.h"
#include "queso/embedding/mesh_operator.h"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/io/io_utilities.h"

#include "queso/tests/cpp/helpers/test_data.hpp"

namespace queso::testing {

    TEST_SUITE("ElementClassifierTestSuite")
    {

        TEST_CASE("ElementClassifierTestSuite - Touch Element Cube Test")
        {


            // Read mesh from STL file
            TriangleMesh triangle_mesh{};
            std::string base_dir = TestDataDirectory();
            IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/primitives/cube_with_cavity.stl");

            // Instatiate brep_operator
            const GeometryTolerance tolerance =
                GeometryTolerance::FromScale({ .length_scale = 4.0, .coordinate_scale = 2.0 });
            embedding::MeshOperator brep_operator(triangle_mesh.View(), tolerance);
            const auto Classify = [&](PointView rLower, PointView rUpper, AabbIntersectionPolicy Policy) {
                const auto candidates = brep_operator.Query().GetAabbCandidates(rLower, rUpper);
                if (brep_operator.Query().IntersectsAabb(candidates, rLower, rUpper, Policy)) {
                    return IntersectionState::trimmed;
                }
                return brep_operator.IsInside(0.5 * (rLower + rUpper)) ? IntersectionState::inside
                                                                       : IntersectionState::outside;
            };

            Vector3d lower_bound = { -2, -2, -2 };
            Vector3d upper_bound = { -1.5, 2, 2 };
            // Touch from outside with tolerance=0.0 is trimmed.
            CHECK((Classify(lower_bound, upper_bound, AabbIntersectionPolicy::Exact)) == (IntersectionState::trimmed));
            // Touch from outside with tolerance>0.0 is outside.
            CHECK(
                (Classify(lower_bound, upper_bound, AabbIntersectionPolicy::SnapEroded)) == (IntersectionState::outside)
            );

            lower_bound = { -1.5, -1.5, -1.5 };
            upper_bound = { -1.4, -1.4, -1.4 };
            // Touch from inside with tolerance=0.0 is trimmed.
            CHECK((Classify(lower_bound, upper_bound, AabbIntersectionPolicy::Exact)) == (IntersectionState::trimmed));
            // Touch from inside with tolerance>0.0 is inside.
            CHECK(
                (Classify(lower_bound, upper_bound, AabbIntersectionPolicy::SnapEroded)) == (IntersectionState::inside)
            );
        }

        TEST_CASE("ElementClassifierTestSuite - Cylinder Element Classifier Test")
        {


            // Read mesh from STL file
            TriangleMesh triangle_mesh{};
            std::string base_dir = TestDataDirectory();
            IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/primitives/cylinder.stl");

            auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
            auto& r_settings = *p_settings;

            auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
            r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
            r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ -1.5, -1.5, -1.0 });
            r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{ 1.5, 1.5, 12.0 });
            r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ -1.0, -1 - 0, 1.0 });
            r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 1.0, 1.0, 1.0 });
            r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
            r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{ 30, 30, 130 });
            r_grid_settings.CheckRequired();

            GridIndexer grid_indexer(r_settings);

            // Get flood fill solution
            embedding::DomainMeshEmbedder embedder(triangle_mesh.View(), grid_indexer);
            const auto classification = embedder.Classify();

            REQUIRE((classification.size()) == (grid_indexer.NumberOfElements()));
            std::ifstream myfile(base_dir + "/embedding/element_classifier_cylinder.txt");
            std::string line;
            for (IndexType i = 0; i < classification.size(); ++i) {
                getline(myfile, line);
                CHECK((static_cast<IndexType>(classification[i])) == (static_cast<IndexType>(std::stoi(line))));
            }
            myfile.close();
        }

        TEST_CASE("ElementClassifierTestSuite - Cube Element Classifier Test")
        {


            // Read mesh from STL file
            TriangleMesh triangle_mesh{};
            std::string base_dir = TestDataDirectory();
            IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/primitives/cube_with_cavity.stl");

            auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
            auto& r_settings = *p_settings;

            auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
            r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
            // Preserve the original 20^3 cells and add one same-sized exterior halo cell on every side.
            r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ -1.65, -1.65, -1.65 });
            r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{ 1.65, 1.65, 1.65 });
            r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ -1.0, -1 - 0, 1.0 });
            r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 1.0, 1.0, 1.0 });
            r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
            r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{ 22, 22, 22 });
            r_grid_settings.CheckRequired();

            GridIndexer grid_indexer(r_settings);

            // Get flood fill solution
            embedding::DomainMeshEmbedder embedder(triangle_mesh.View(), grid_indexer);
            const auto classification = embedder.Classify();

            REQUIRE((classification.size()) == (grid_indexer.NumberOfElements()));

            std::ifstream myfile(base_dir + "/embedding/element_classifier_cube.txt");
            std::string line;
            constexpr Vector3i reference_number_of_elements{ 20, 20, 20 };
            IndexType reference_index = 0;
            for (IndexType k = 0; k < reference_number_of_elements[2]; ++k) {
                for (IndexType j = 0; j < reference_number_of_elements[1]; ++j) {
                    for (IndexType i = 0; i < reference_number_of_elements[0]; ++i) {
                        getline(myfile, line);
                        const IndexType index = grid_indexer.GetVectorIndexFromMatrixIndices(i + 1, j + 1, k + 1);
                        CHECK(
                            (static_cast<IndexType>(classification[index])) == (static_cast<IndexType>(std::stoi(line)))
                        );
                        ++reference_index;
                    }
                }
            }
            CHECK((reference_index) == (20UL * 20UL * 20UL));
            myfile.close();
        }

        TEST_CASE("ElementClassifierTestSuite - Elephant Element Classifier Test")
        {


            // Read mesh from STL file
            TriangleMesh triangle_mesh{};
            std::string base_dir = TestDataDirectory();
            IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/primitives/elephant.stl");

            auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
            auto& r_settings = *p_settings;

            auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
            r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
            r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ -0.4, -0.6, -0.35 });
            r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{ 0.4, 0.6, 0.35 });
            r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ -1.0, -1 - 0, 1.0 });
            r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 1.0, 1.0, 1.0 });
            r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
            r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{ 16, 24, 14 });
            r_grid_settings.CheckRequired();

            GridIndexer grid_indexer(r_settings);

            // Get flood fill solution
            embedding::DomainMeshEmbedder embedder(triangle_mesh.View(), grid_indexer);
            const auto classification = embedder.Classify();

            REQUIRE((classification.size()) == (5376));
            std::ifstream myfile(base_dir + "/embedding/element_classifier_elephant.txt");
            std::string line;
            for (IndexType i = 0; i < classification.size(); ++i) {
                getline(myfile, line);
                CHECK((static_cast<IndexType>(classification[i])) == (static_cast<IndexType>(std::stoi(line))));
            }
            myfile.close();
        }

        TEST_CASE("ElementClassifierTestSuite - Bunny Element Classifier Test")
        {


            // Read mesh from STL file
            TriangleMesh triangle_mesh{};
            std::string base_dir = TestDataDirectory();
            IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/primitives/stanford_bunny.stl");

            auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
            auto& r_settings = *p_settings;

            auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
            r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
            r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ -24, -43, 5 });
            r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{ 85, 46, 115 });
            r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ -1.0, -1 - 0, 1.0 });
            r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 1.0, 1.0, 1.0 });
            r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
            r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{ 36, 30, 40 });
            r_grid_settings.CheckRequired();

            GridIndexer grid_indexer(r_settings);

            // Get flood fill solution
            embedding::DomainMeshEmbedder embedder(triangle_mesh.View(), grid_indexer);
            const auto classification = embedder.Classify();

            REQUIRE((classification.size()) == (grid_indexer.NumberOfElements()));
            std::ifstream myfile(base_dir + "/embedding/element_classifier_bunny.txt");
            std::string line;
            for (IndexType i = 0; i < classification.size(); ++i) {
                getline(myfile, line);
                CHECK((static_cast<IndexType>(classification[i])) == (static_cast<IndexType>(std::stoi(line))));
            }
            myfile.close();
        }

    }  // TEST_SUITE

} // namespace queso::testing
