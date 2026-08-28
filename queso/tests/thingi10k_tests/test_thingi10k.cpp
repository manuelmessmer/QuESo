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
#ifdef _OPENMP
#include <omp.h>
#endif

//// STL includes
#include <cstdint>
#include <fstream>
#include <string>

//// Project includes
#include "queso/containers/background_grid.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/containers/untrimmed_element.hpp"
#include "queso/embedding/boundary_mesh_embedder.h"
#include "queso/embedding/domain_mesh_embedder.h"
#include "queso/embedding/flood_fill.h"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/io/io_utilities.h"
#include "queso/utilities/mesh_utilities.h"

namespace queso::Testing {
namespace {

    using BackgroundGridType = BackgroundGrid<IntegrationPoint, BoundaryIntegrationPoint>;
    using UntrimmedElementType = BackgroundGridType::UntrimmedElementType;

    struct ActiveCellBuilder
    {
        static constexpr BackgroundGridType::ElementFilter Builds = BackgroundGridType::ElementFilter::untrimmed;

        [[nodiscard]] std::optional<UntrimmedElementType> Build(IndexType CellIndex, const ElementBounds& rBounds)
        { return UntrimmedElementType(CellIndex + 1, rBounds); }
    };

    [[nodiscard]] Unique<MainDictionaryType>
        MakeSettings(const std::string& rFilename, const BoundingBoxType& rBounds, const Vector3i& rNumberOfElements)
    {
        auto p_settings = factories::CreateSettings();
        p_settings->operator[](MainSettings::general_settings).SetValue(GeneralSettings::input_filename, rFilename);
        auto& r_grid_settings = (*p_settings)[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, rBounds.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, rBounds.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, rBounds.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, rBounds.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, rNumberOfElements);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        p_settings->CheckRequired();
        return p_settings;
    }

    [[nodiscard]] std::vector<std::string> GetFilenames(
        int ArgumentCount,
        char** pArguments,
        IndexType NumberOfElementsMin,
        IndexType NumberOfElementsMax,
        std::string_view TestName
    )
    {
        QuESo_ERROR_IF(ArgumentCount != 6)
            << "Please provide following arguments: -- single/small_set/large_set Filename/Directory n_min n_max\n";

        std::vector<std::string> filenames;
        filenames.reserve(5000);
        const std::string option = pArguments[1];
        if (option == "single") {
            filenames.emplace_back(pArguments[2]);
            QuESo_INFO << TestName << " :: Testing single STL: " << filenames.front()
                       << " with n_min: " << NumberOfElementsMin << ", n_max: " << NumberOfElementsMax << ".\n";
            return filenames;
        }

        QuESo_ERROR_IF(option != "small_set" && option != "large_set")
            << "First argument must be: 'single', 'small_set', or 'large_set'. Provided: " << option;
        const std::string list_filename = std::string(pArguments[5]) + "/model_ids_" + option + ".txt";
        std::ifstream file(list_filename);
        QuESo_ERROR_IF(!file) << "Could not open model list: " << list_filename << '\n';
        std::string model_id;
        while (std::getline(file, model_id)) { filenames.emplace_back(std::string(pArguments[2]) + model_id + ".stl"); }
        QuESo_INFO << "Testing '" << option << "' containing " << filenames.size()
                   << " STLs with n_min: " << NumberOfElementsMin << ", n_max: " << NumberOfElementsMax << ".\n";
        return filenames;
    }

    struct GridConfiguration
    {
        BoundingBoxType bounds;
        Vector3i number_of_elements;
    };

    [[nodiscard]] GridConfiguration MakeGridConfiguration(
        const TriangleMeshView& rMesh,
        IndexType NumberOfElementsMin,
        IndexType NumberOfElementsMax
    )
    {
        const auto [lower, upper] = MeshUtilities::BoundingBox(rMesh);
        const PointType delta = upper - lower;
        const double cell_size =
            1.2 * std::min(Math::Max(delta) / NumberOfElementsMax, Math::Min(delta) / NumberOfElementsMin);
        const PointType grid_lower = lower - 0.1 * delta;
        Vector3i number_of_elements{};
        PointType grid_upper{};
        for (IndexType axis = 0; axis < 3; ++axis) {
            number_of_elements[axis] = static_cast<IndexType>(std::ceil(1.2 * delta[axis] / cell_size));
            grid_upper[axis] = grid_lower[axis] + number_of_elements[axis] * cell_size;
        }
        return { { grid_lower, grid_upper }, number_of_elements };
    }

    [[nodiscard]] BackgroundGridType
        MakeActiveGrid(const MainDictionaryType& rSettings, const embedding::ElementStates& rStates)
    {
        BackgroundGridType grid(rSettings);
        ActiveCellBuilder builder;
        grid.ReserveElements(
            std::ranges::count_if(rStates, [](IntersectionState State) { return State != IntersectionState::outside; }),
            0
        );
        for (IndexType cell_index = 0; cell_index < rStates.size(); ++cell_index) {
            if (rStates[cell_index] == IntersectionState::outside) { continue; }
            QuESo_ERROR_IF(!grid.MakeElement(builder, cell_index))
                << "Failed to construct active placeholder for cell " << cell_index << ".\n";
        }
        grid.LockElements();
        return grid;
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(Thingi10KTestSuite)

BOOST_AUTO_TEST_CASE(STLEmbeddingTest)
{
    QuESo_INFO << "Testing :: Test Thingi10k :: STL embedding test.\n";
    Timer timer;
    const auto& r_master_suite = boost::unit_test::framework::master_test_suite();
    const IndexType number_of_elements_min = static_cast<IndexType>(std::stoi(r_master_suite.argv[3]));
    const IndexType number_of_elements_max = static_cast<IndexType>(std::stoi(r_master_suite.argv[4]));
    const auto filenames = GetFilenames(
        r_master_suite.argc,
        r_master_suite.argv,
        number_of_elements_min,
        number_of_elements_max,
        "Thingi10KSTLEmbeddingTest"
    );

    std::int64_t successful_tests = 0;
    double max_relative_volume_error = 0.0;
    double max_relative_area_error = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : successful_tests) schedule(dynamic)
#endif
    for (std::int64_t raw_file_index = 0; raw_file_index < static_cast<std::int64_t>(filenames.size());
         ++raw_file_index) {
        const std::string& r_filename = filenames[static_cast<IndexType>(raw_file_index)];
        TriangleMesh mesh;
        IO::ReadMeshFromSTL(mesh, r_filename);

        const GridConfiguration configuration =
            MakeGridConfiguration(mesh.View(), number_of_elements_min, number_of_elements_max);
        auto p_settings = MakeSettings(r_filename, configuration.bounds, configuration.number_of_elements);
        BackgroundGridType active_grid(*p_settings);
        const GridIndexer& r_grid_indexer = active_grid.GetGridIndexer();
        embedding::DomainMeshEmbedder domain_embedder(mesh.View(), r_grid_indexer);
        const auto states = domain_embedder.Classify();
        active_grid = MakeActiveGrid(*p_settings, states);

        double represented_volume = 0.0;
        constexpr IndexType minimum_boundary_triangles = 10;
        for (IndexType cell_index = 0; cell_index < states.size(); ++cell_index) {
            const BoundingBoxType bounds = r_grid_indexer.GetBoundingBoxXYZFromIndex(cell_index);
            if (states[cell_index] == IntersectionState::trimmed) {
                represented_volume += MeshUtilities::Volume(
                    domain_embedder.MakeTrimmedDomain(cell_index, minimum_boundary_triangles).GetBoundaryMesh()
                );
            } else if (states[cell_index] == IntersectionState::inside) {
                const PointType delta = bounds.upper - bounds.lower;
                represented_volume += delta[0] * delta[1] * delta[2];
            }
        }

        embedding::BoundaryMeshEmbedder boundary_embedder(mesh.View(), active_grid);
        double represented_area = 0.0;
        IndexType consumed_sections = 0;
        for (const IndexType cell_index : boundary_embedder.GetParentCellIndices()) {
            auto section = boundary_embedder.TakeBoundarySection(cell_index);
            represented_area += MeshUtilities::Area(section.View());
            ++consumed_sections;
        }
        const double reference_area = MeshUtilities::Area(mesh.View());
        const double area_error = std::abs(represented_area - reference_area) / reference_area;
        const double reference_volume = MeshUtilities::Volume(mesh.View());
        const double volume_error = std::abs(represented_volume - reference_volume) / reference_volume;
#ifdef _OPENMP
#pragma omp critical(thingi10k_error_maximum)
#endif
        {
            max_relative_volume_error = std::max(max_relative_volume_error, volume_error);
            max_relative_area_error = std::max(max_relative_area_error, area_error);
        }

        constexpr double volume_tolerance = 1e-8;
        constexpr double area_tolerance = 1e-10;
        BOOST_CHECK_LT(volume_error, volume_tolerance);
        BOOST_CHECK_LT(area_error, area_tolerance);
        if (volume_error < volume_tolerance && area_error < area_tolerance) {
            ++successful_tests;
        } else {
            QuESo_INFO << "Test failed for filename: " << r_filename << ", volume error: " << volume_error
                       << ", area error: " << area_error
                       << ", listed sections: " << boundary_embedder.GetParentCellIndices().size()
                       << ", consumed sections: " << consumed_sections << '\n';
        }
    }

    QuESo_INFO << "Successful tests: " << successful_tests << '\n';
    QuESo_INFO << "Maximum relative volume error: " << max_relative_volume_error << '\n';
    QuESo_INFO << "Maximum relative area error: " << max_relative_area_error << '\n';
    QuESo_INFO << "Elapsed time: " << timer.Measure() << " sec\n";
}

BOOST_AUTO_TEST_CASE(ElementClassificationTest)
{
    QuESo_INFO << "Testing :: Test Thingi10k :: Element classification.\n";
    Timer timer;
    const auto& r_master_suite = boost::unit_test::framework::master_test_suite();
    const IndexType number_of_elements_min = static_cast<IndexType>(std::stoi(r_master_suite.argv[3]));
    const IndexType number_of_elements_max = static_cast<IndexType>(std::stoi(r_master_suite.argv[4]));
    const auto filenames = GetFilenames(
        r_master_suite.argc,
        r_master_suite.argv,
        number_of_elements_min,
        number_of_elements_max,
        "Thingi10KElementClassificationTest"
    );

    IndexType successful_tests = 0;
    for (const std::string& r_filename : filenames) {
        TriangleMesh mesh;
        IO::ReadMeshFromSTL(mesh, r_filename);
        const GridConfiguration configuration =
            MakeGridConfiguration(mesh.View(), number_of_elements_min, number_of_elements_max);
        const auto p_settings = MakeSettings(r_filename, configuration.bounds, configuration.number_of_elements);
        GridIndexer grid_indexer(*p_settings);
        embedding::DomainMeshEmbedder domain_embedder(mesh.View(), grid_indexer);
        const auto states = domain_embedder.Classify();
        embedding::DomainMeshEmbedder serial_domain_embedder(mesh.View(), grid_indexer);
        const auto serial_states = serial_domain_embedder.Classify({ .partition_count = 1 });

        QuESo_ERROR_IF(states.size() != serial_states.size())
            << "Classification size mismatch for " << r_filename << ".\n";
        for (IndexType cell_index = 0; cell_index < states.size(); ++cell_index) {
            QuESo_ERROR_IF(states[cell_index] != serial_states[cell_index])
                << "Partition-dependent classification for " << r_filename << " at cell " << cell_index << ".\n";
        }
        ++successful_tests;
    }

    QuESo_INFO << "Successful tests: " << successful_tests << '\n';
    QuESo_INFO << "Elapsed time: " << timer.Measure() << " sec\n";
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing

