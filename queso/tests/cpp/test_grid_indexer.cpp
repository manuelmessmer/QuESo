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
#include <array>
#include <initializer_list>
#include <limits>
#include <set>
#include <vector>
//// Project includes
#include "queso/containers/grid_indexer.hpp"
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/utilities/math_utilities.hpp"

namespace queso {
namespace Testing {

    namespace {

        [[nodiscard]] Unique<Dictionary<queso::key::MainValuesTypeTag>>
            MakeGridSettings(const BoundingBoxType& rBounds, const Vector3i& rNumberOfElements)
        {
            auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
            auto& r_grid_settings = (*p_settings)[MainSettings::background_grid_settings];
            r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
            r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, rBounds.lower);
            r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, rBounds.upper);
            r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, rBounds.lower);
            r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, rBounds.upper);
            r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
            r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, rNumberOfElements);
            r_grid_settings.CheckRequired();
            return p_settings;
        }

    }  // namespace

    BOOST_AUTO_TEST_SUITE(GridIndexerTestSuite)

    BOOST_AUTO_TEST_CASE(GridIndexerBSplineMeshTest)
    {
        QuESo_INFO << "Testing :: Test Grid Indexer :: BSpline Mesh" << std::endl;

        const BoundingBoxType bounds_xyz = MakeBox({ -1.0, -0.5, 1.0 }, { 5.0, 10.5, 13.0 });
        const BoundingBoxType bounds_uvw = MakeBox({ -1.0, -0.5, 1.0 }, { 5.0, 10.5, 13.0 });

        const Vector3i number_of_elements{ 5, 10, 7 };

        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_settings = *p_settings;

        auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);
        r_grid_settings.CheckRequired();

        GridIndexer grid_indexer(r_settings);
        const auto delta = (bounds_xyz.upper - bounds_xyz.lower);
        double volume = 0.0;
        for (IndexType i = 0; i < number_of_elements[0]; ++i) {
            for (IndexType j = 0; j < number_of_elements[1]; ++j) {
                for (IndexType k = 0; k < number_of_elements[2]; ++k) {
                    IndexType index = grid_indexer.GetVectorIndexFromMatrixIndices(i, j, k);
                    auto indices = grid_indexer.GetMatrixIndicesFromVectorIndex(index);
                    QuESo_CHECK_EQUAL(i, indices[0]);
                    QuESo_CHECK_EQUAL(j, indices[1]);
                    QuESo_CHECK_EQUAL(k, indices[2]);

                    auto box_1_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(index);
                    auto box_2_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(i, j, k);
                    auto box_3_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(indices);

                    QuESo_CHECK_POINT_NEAR(box_1_xyz.lower, box_2_xyz.lower, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_1_xyz.lower, box_3_xyz.lower, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_1_xyz.upper, box_2_xyz.upper, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_1_xyz.upper, box_3_xyz.upper, 1e-12);

                    auto box_1_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(index);
                    auto box_2_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(i, j, k);
                    auto box_3_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(indices);

                    QuESo_CHECK_POINT_NEAR(box_1_uvw.lower, box_1_xyz.lower, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_1_uvw.upper, box_1_xyz.upper, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_2_uvw.lower, box_2_xyz.lower, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_2_uvw.upper, box_2_xyz.upper, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_3_uvw.lower, box_3_xyz.lower, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_3_uvw.upper, box_3_xyz.upper, 1e-12);

                    auto delta_box = (box_1_xyz.upper - box_1_xyz.lower);
                    QuESo_CHECK_LT(std::abs(delta[0] / number_of_elements[0] - delta_box[0]) / delta_box[0], 1e-12);
                    QuESo_CHECK_LT(std::abs(delta[1] / number_of_elements[1] - delta_box[1]) / delta_box[1], 1e-12);
                    QuESo_CHECK_LT(std::abs(delta[2] / number_of_elements[2] - delta_box[2]) / delta_box[2], 1e-12);

                    volume += delta_box[0] * delta_box[1] * delta_box[2];

                    IndexType index_1 = grid_indexer.GetVectorIndexFromMatrixIndices(i, j, k);
                    IndexType index_2 = grid_indexer.GetVectorIndexFromMatrixIndices(indices);
                    QuESo_CHECK_EQUAL(index_1, index);
                    QuESo_CHECK_EQUAL(index_2, index);
                }
            }
        }

        double volume_ref = delta[0] * delta[1] * delta[2];
        QuESo_CHECK_LT(std::abs(volume - volume_ref) / volume_ref, 1e-12);
    }


    BOOST_AUTO_TEST_CASE(GridIndexerFEMeshTest)
    {
        QuESo_INFO << "Testing :: Test Grid Indexer :: Hexahedral FE Mesh" << std::endl;

        const BoundingBoxType bounds_xyz = MakeBox({ -1.0, -0.5, 1.0 }, { 5.0, 10.5, 13.0 });
        const BoundingBoxType bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });

        const Vector3i number_of_elements{ 5, 10, 7 };

        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_settings = *p_settings;

        auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);
        r_grid_settings.CheckRequired();

        GridIndexer grid_indexer(r_settings);
        const auto delta = (bounds_xyz.upper - bounds_xyz.lower);
        double volume = 0.0;
        for (IndexType i = 0; i < number_of_elements[0]; ++i) {
            for (IndexType j = 0; j < number_of_elements[1]; ++j) {
                for (IndexType k = 0; k < number_of_elements[2]; ++k) {
                    IndexType index = grid_indexer.GetVectorIndexFromMatrixIndices(i, j, k);
                    auto indices = grid_indexer.GetMatrixIndicesFromVectorIndex(index);
                    QuESo_CHECK_EQUAL(i, indices[0]);
                    QuESo_CHECK_EQUAL(j, indices[1]);
                    QuESo_CHECK_EQUAL(k, indices[2]);

                    auto box_1_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(index);
                    auto box_2_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(i, j, k);
                    auto box_3_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(indices);

                    QuESo_CHECK_POINT_NEAR(box_1_xyz.lower, box_2_xyz.lower, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_1_xyz.lower, box_3_xyz.lower, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_1_xyz.upper, box_2_xyz.upper, 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_1_xyz.upper, box_3_xyz.upper, 1e-12);

                    auto box_1_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(index);
                    auto box_2_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(i, j, k);
                    auto box_3_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(indices);

                    QuESo_CHECK_POINT_NEAR(box_1_uvw.lower, PointType({ -1.0, -1.0, -1.0 }), 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_1_uvw.upper, PointType({ 1.0, 1.0, 1.0 }), 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_2_uvw.lower, PointType({ -1.0, -1.0, -1.0 }), 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_2_uvw.upper, PointType({ 1.0, 1.0, 1.0 }), 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_3_uvw.lower, PointType({ -1.0, -1.0, -1.0 }), 1e-12);
                    QuESo_CHECK_POINT_NEAR(box_3_uvw.upper, PointType({ 1.0, 1.0, 1.0 }), 1e-12);

                    auto delta_box = (box_1_xyz.upper - box_1_xyz.lower);
                    QuESo_CHECK_LT(std::abs(delta[0] / number_of_elements[0] - delta_box[0]) / delta_box[0], 1e-12);
                    QuESo_CHECK_LT(std::abs(delta[1] / number_of_elements[1] - delta_box[1]) / delta_box[1], 1e-12);
                    QuESo_CHECK_LT(std::abs(delta[2] / number_of_elements[2] - delta_box[2]) / delta_box[2], 1e-12);

                    volume += delta_box[0] * delta_box[1] * delta_box[2];

                    IndexType index_1 = grid_indexer.GetVectorIndexFromMatrixIndices(i, j, k);
                    IndexType index_2 = grid_indexer.GetVectorIndexFromMatrixIndices(indices);
                    QuESo_CHECK_EQUAL(index_1, index);
                    QuESo_CHECK_EQUAL(index_2, index);
                }
            }
        }

        double volume_ref = delta[0] * delta[1] * delta[2];
        QuESo_CHECK_LT(std::abs(volume - volume_ref) / volume_ref, 1e-12);
    }

    bool contains(std::vector<IndexType>& v, IndexType test_value)
    {
        if (std::find(v.begin(), v.end(), test_value) != v.end()) { return true; }
        return false;
    }

    BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingGlobalXTest)
    {
        QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Global Partition X" << std::endl;

        const BoundingBoxType bounds_xyz = MakeBox({ -1.0, -0.5, 1.0 }, { 5.0, 10.5, 13.0 });
        const BoundingBoxType bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });

        const Vector3i number_of_elements{ 3, 4, 5 };

        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_settings = *p_settings;

        auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);
        r_grid_settings.CheckRequired();

        GridIndexer grid_indexer(r_settings);

        // Walk foward
        IndexType index = 0;
        GridIndexer::IndexInfo index_info{};
        IndexType i = 0;
        while (index < grid_indexer.NumberOfElements() - 1) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if ((index + 1) % 3 == 0) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::x_forward>(index);
            ++i;
            QuESo_CHECK_EQUAL(index, i);
            QuESo_CHECK_EQUAL(index_info, index_info_ref);
        }
        QuESo_CHECK_EQUAL(index, 59);
        std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::x_forward>(index);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 59);
        --i;

        // Walk backwards
        while (index > 0) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if ((index) % 3 == 0) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::x_backward>(index);
            QuESo_CHECK_EQUAL(index, i);
            QuESo_CHECK_EQUAL(index_info, index_info_ref);
            --i;
        }
        QuESo_CHECK_EQUAL(index, 0);
        std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::x_backward>(index);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 0);
    }

    BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingGlobalYTest)
    {
        QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Global Partition Y" << std::endl;

        const BoundingBoxType bounds_xyz = MakeBox({ -1.0, -0.5, 1.0 }, { 5.0, 10.5, 13.0 });
        const BoundingBoxType bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });

        const Vector3i number_of_elements{ 3, 4, 5 };

        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_settings = *p_settings;

        auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);
        r_grid_settings.CheckRequired();

        GridIndexer grid_indexer(r_settings);
        // Walk foward
        IndexType index = 0;
        GridIndexer::IndexInfo index_info;
        IndexType i = 1;
        std::vector<IndexType> order_y = { 0,  3,  6,  9,  1,  4,  7,  10, 2,  5,  8,  11, 12, 15, 18,
                                           21, 13, 16, 19, 22, 14, 17, 20, 23, 24, 27, 30, 33, 25, 28,
                                           31, 34, 26, 29, 32, 35, 36, 39, 42, 45, 37, 40, 43, 46, 38,
                                           41, 44, 47, 48, 51, 54, 57, 49, 52, 55, 58, 50, 53, 56, 59 };
        std::vector<IndexType> local_ends_y = { 9, 10, 11, 21, 22, 23, 33, 34, 35, 45, 46, 47, 57, 58 };
        while (index < grid_indexer.NumberOfElements() - 1) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends_y, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::y_forward>(index);
            QuESo_CHECK_EQUAL(order_y[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
            ++i;
        }

        QuESo_CHECK_EQUAL(index, 59);
        std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::y_forward>(index);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 59);
        --i;

        local_ends_y = { 0, 1, 2, 12, 13, 14, 24, 25, 26, 36, 37, 38, 48, 49, 50 };
        while (index > 0) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends_y, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::y_backward>(index);
            --i;
            QuESo_CHECK_EQUAL(order_y[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
        }
        QuESo_CHECK_EQUAL(index, 0);
        std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::z_backward>(index);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 0);
    }


    BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingGlobalZTest)
    {
        QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Global Partition Z" << std::endl;

        const BoundingBoxType bounds_xyz = MakeBox({ -1.0, -0.5, 1.0 }, { 5.0, 10.5, 13.0 });
        const BoundingBoxType bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });

        const Vector3i number_of_elements{ 3, 4, 5 };

        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_settings = *p_settings;

        auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);
        r_grid_settings.CheckRequired();

        GridIndexer grid_indexer(r_settings);
        // Walk foward
        IndexType index = 0;
        GridIndexer::IndexInfo index_info;
        IndexType i = 1;
        std::vector<IndexType> order_z = {
            0, 12, 24, 36, 48, 1, 13, 25, 37, 49, 2,  14, 26, 38, 50, 3,  15, 27, 39, 51,
            4, 16, 28, 40, 52, 5, 17, 29, 41, 53, 6,  18, 30, 42, 54, 7,  19, 31, 43, 55,
            8, 20, 32, 44, 56, 9, 21, 33, 45, 57, 10, 22, 34, 46, 58, 11, 23, 35, 47, 59
        };

        std::vector<IndexType> local_ends_z = { 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58 };
        while (index < grid_indexer.NumberOfElements() - 1) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends_z, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::z_forward>(index);
            QuESo_CHECK_EQUAL(order_z[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
            ++i;
        }

        QuESo_CHECK_EQUAL(index, 59);
        std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::z_forward>(index);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 59);
        --i;

        // Walk backwards
        local_ends_z = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
        while (index > 0) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends_z, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::z_backward>(index);
            --i;
            QuESo_CHECK_EQUAL(order_z[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
        }
        QuESo_CHECK_EQUAL(index, 0);
        std::tie(index, index_info) = grid_indexer.GetNextIndex<GridIndexer::Direction::z_backward>(index);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 0);
    }

    BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingLocalXTest)
    {
        QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Local Partition X" << std::endl;

        const BoundingBoxType bounds_xyz = MakeBox({ -1.0, -0.5, 1.0 }, { 5.0, 10.5, 13.0 });
        const BoundingBoxType bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });

        const Vector3i number_of_elements{ 3, 4, 5 };

        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_settings = *p_settings;

        auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);
        r_grid_settings.CheckRequired();

        GridIndexer grid_indexer(r_settings);

        PartitionBoxType local_partition_1 = std::make_pair(Vector3i({ 1, 1, 1 }), Vector3i({ 2, 3, 3 }));
        // Walk foward
        IndexType index = 16;
        GridIndexer::IndexInfo index_info;
        IndexType i = 1;
        std::vector<IndexType> order = { 16, 17, 19, 20, 22, 23, 28, 29, 31, 32, 34, 35, 40, 41, 43, 44, 46, 47 };
        std::vector<IndexType> local_ends = { 17, 20, 23, 29, 32, 35, 41, 44, 47 };
        while (index < 47) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::x_forward>(index, local_partition_1);
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
            ++i;
        }

        QuESo_CHECK_EQUAL(index, 47);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::x_forward>(index, local_partition_1);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 47);
        --i;

        // Walk backwards
        local_ends = { 16, 19, 22, 28, 31, 34, 40, 43, 46 };
        while (index > 16) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::x_backward>(index, local_partition_1);
            --i;
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
        }
        QuESo_CHECK_EQUAL(index, 16);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::x_backward>(index, local_partition_1);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 16);

        PartitionBoxType local_partition_2 = std::make_pair(Vector3i({ 0, 0, 0 }), Vector3i({ 1, 1, 2 }));
        // Walk foward
        index = 0;
        i = 1;
        order = { 0, 1, 3, 4, 12, 13, 15, 16, 24, 25, 27, 28 };
        local_ends = { 1, 4, 13, 16, 25, 28 };
        while (index < 28) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::x_forward>(index, local_partition_2);
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
            ++i;
        }

        QuESo_CHECK_EQUAL(index, 28);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::x_forward>(index, local_partition_2);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 28);
        --i;

        // Walk backwards
        local_ends = { 0, 3, 12, 15, 24, 27 };
        while (index > 0) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::x_backward>(index, local_partition_2);
            --i;
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
        }
        QuESo_CHECK_EQUAL(index, 0);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::x_backward>(index, local_partition_2);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 0);
    }

    BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingLocalYTest)
    {
        QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Local Partition Y" << std::endl;

        const BoundingBoxType bounds_xyz = MakeBox({ -1.0, -0.5, 1.0 }, { 5.0, 10.5, 13.0 });
        const BoundingBoxType bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });

        const Vector3i number_of_elements{ 3, 4, 5 };

        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_settings = *p_settings;

        auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);
        r_grid_settings.CheckRequired();

        GridIndexer grid_indexer(r_settings);

        PartitionBoxType local_partition_1 = std::make_pair(Vector3i({ 1, 1, 1 }), Vector3i({ 2, 3, 3 }));
        // Walk foward
        IndexType index = 16;
        GridIndexer::IndexInfo index_info;
        IndexType i = 1;
        std::vector<IndexType> order = { 16, 19, 22, 17, 20, 23, 28, 31, 34, 29, 32, 35, 40, 43, 46, 41, 44, 47 };
        std::vector<IndexType> local_ends = { 22, 23, 34, 35, 46, 47 };
        while (index < 47) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::y_forward>(index, local_partition_1);
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
            ++i;
        }

        QuESo_CHECK_EQUAL(index, 47);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::y_forward>(index, local_partition_1);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 47);
        --i;

        // Walk backwards
        local_ends = { 16, 17, 28, 29, 40, 41 };
        while (index > 16) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::y_backward>(index, local_partition_1);
            --i;
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
        }
        QuESo_CHECK_EQUAL(index, 16);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::y_backward>(index, local_partition_1);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 16);


        PartitionBoxType local_partition_2 = std::make_pair(Vector3i({ 0, 0, 0 }), Vector3i({ 1, 1, 2 }));
        // Walk foward
        index = 0;
        i = 1;
        order = { 0, 3, 1, 4, 12, 15, 13, 16, 24, 27, 25, 28 };
        local_ends = { 3, 4, 15, 16, 27, 28 };
        while (index < 28) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::y_forward>(index, local_partition_2);
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
            ++i;
        }

        QuESo_CHECK_EQUAL(index, 28);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::y_forward>(index, local_partition_2);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 28);
        --i;

        // Walk backwards
        local_ends = { 0, 1, 12, 13, 24, 25 };
        while (index > 0) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::y_backward>(index, local_partition_2);
            --i;
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
        }
        QuESo_CHECK_EQUAL(index, 0);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::y_backward>(index, local_partition_2);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 0);
    }


    BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingLocalZTest)
    {
        QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Local Partition Z" << std::endl;

        const BoundingBoxType bounds_xyz = MakeBox({ -1.0, -0.5, 1.0 }, { 5.0, 10.5, 13.0 });
        const BoundingBoxType bounds_uvw = MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });

        const Vector3i number_of_elements{ 3, 4, 5 };

        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_settings = *p_settings;

        auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);
        r_grid_settings.CheckRequired();

        GridIndexer grid_indexer(r_settings);

        PartitionBoxType local_partition_1 = std::make_pair(Vector3i({ 1, 1, 1 }), Vector3i({ 2, 3, 3 }));
        // Walk foward
        IndexType index = 16;
        GridIndexer::IndexInfo index_info;
        IndexType i = 1;
        std::vector<IndexType> order = { 16, 28, 40, 17, 29, 41, 19, 31, 43, 20, 32, 44, 22, 34, 46, 23, 35, 47 };
        std::vector<IndexType> local_ends = { 40, 41, 43, 44, 46, 47 };
        while (index < 47) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::z_forward>(index, local_partition_1);
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
            ++i;
        }

        QuESo_CHECK_EQUAL(index, 47);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::z_forward>(index, local_partition_1);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 47);
        --i;

        // Walk backwards
        local_ends = { 16, 17, 19, 20, 22, 23 };
        while (index > 16) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::z_backward>(index, local_partition_1);
            --i;
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
        }
        QuESo_CHECK_EQUAL(index, 16);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::z_backward>(index, local_partition_1);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 16);

        PartitionBoxType local_partition_2 = std::make_pair(Vector3i({ 0, 0, 0 }), Vector3i({ 1, 1, 2 }));
        // Walk foward
        index = 0;
        i = 1;
        order = { 0, 12, 24, 1, 13, 25, 3, 15, 27, 4, 16, 28 };
        local_ends = { 27, 28, 24, 25 };
        while (index < 28) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::z_forward>(index, local_partition_2);
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
            ++i;
        }

        QuESo_CHECK_EQUAL(index, 28);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::z_forward>(index, local_partition_2);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 28);
        --i;

        // Walk backwards
        local_ends = { 0, 1, 3, 4 };
        while (index > 0) {
            GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
            if (contains(local_ends, index)) { index_info_ref = GridIndexer::IndexInfo::local_end; }
            std::tie(index, index_info) =
                grid_indexer.GetNextIndex<GridIndexer::Direction::z_backward>(index, local_partition_2);
            --i;
            QuESo_CHECK_EQUAL(order[i], index) QuESo_CHECK_EQUAL(index_info, index_info_ref);
        }
        QuESo_CHECK_EQUAL(index, 0);
        std::tie(index, index_info) =
            grid_indexer.GetNextIndex<GridIndexer::Direction::z_backward>(index, local_partition_2);
        QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
        QuESo_CHECK_EQUAL(index, 0);
    }


    BOOST_AUTO_TEST_CASE(RejectsCellsBelowSnapToleranceScale)
    {
        auto p_settings = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 1e-15, 1e-15, 1e-15 }), { 1, 1, 1 });
        BOOST_CHECK_THROW((void)GridIndexer(*p_settings), queso::Exception);
    }

    BOOST_AUTO_TEST_CASE(GeometryToleranceDependsOnCellScaleNotCellCount)
    {
        auto one_cell_settings = MakeGridSettings(MakeBox({ -0.5, -0.5, -0.5 }, { 0.5, 0.5, 0.5 }), { 1, 1, 1 });
        auto two_cell_settings = MakeGridSettings(MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 }), { 2, 2, 2 });
        const GridIndexer one_cell(*one_cell_settings);
        const GridIndexer two_cells(*two_cell_settings);

        QuESo_CHECK_EQUAL(
            one_cell.GetGeometryTolerance().SnapDistance(), two_cells.GetGeometryTolerance().SnapDistance()
        );
        QuESo_CHECK_EQUAL(one_cell.GetGeometryTolerance().ZeroLength(), two_cells.GetGeometryTolerance().ZeroLength());
    }

    BOOST_AUTO_TEST_CASE(GeometryToleranceAccountsForCoordinateMagnitude)
    {
        auto origin_settings = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }), { 1, 1, 1 });
        auto translated_settings =
            MakeGridSettings(MakeBox({ 1e12, 1e12, 1e12 }, { 1e12 + 1.0, 1e12 + 1.0, 1e12 + 1.0 }), { 1, 1, 1 });
        const GridIndexer origin(*origin_settings);
        const GridIndexer translated(*translated_settings);
        const double expected = 32.0 * std::numeric_limits<double>::epsilon() * (1e12 + 1.0);

        QuESo_CHECK_EQUAL(translated.GetGeometryTolerance().SnapDistance(), expected);
        QuESo_CHECK_GT(translated.GetGeometryTolerance().SnapDistance(), origin.GetGeometryTolerance().SnapDistance());
    }

    BOOST_AUTO_TEST_CASE(GeometryToleranceUsesLargestCellDimension)
    {
        auto p_settings = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 4.0, 2.0, 1.0 }), { 1, 1, 1 });
        const GridIndexer grid_indexer(*p_settings);

        const GeometryTolerance unit = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
        const auto& r_tolerance = grid_indexer.GetGeometryTolerance();
        QuESo_CHECK_EQUAL(r_tolerance.SnapDistance(), 4.0 * unit.SnapDistance());
        QuESo_CHECK_EQUAL(r_tolerance.ZeroLength(), 4.0 * unit.ZeroLength());
        QuESo_CHECK_EQUAL(r_tolerance.ZeroArea(), 16.0 * unit.ZeroArea());
        QuESo_CHECK_EQUAL(r_tolerance.ZeroVolume(), 64.0 * unit.ZeroVolume());
    }

    BOOST_AUTO_TEST_CASE(RejectsMalformedGridConfiguration)
    {
        auto zero_count = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }), { 0, 1, 1 });
        BOOST_CHECK_THROW((void)GridIndexer(*zero_count), queso::Exception);

        auto inverted = MakeGridSettings(MakeBox({ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 1.0 }), { 1, 1, 1 });
        BOOST_CHECK_THROW((void)GridIndexer(*inverted), queso::Exception);

        auto non_finite = MakeGridSettings(
            MakeBox({ 0.0, 0.0, 0.0 }, { std::numeric_limits<double>::infinity(), 1.0, 1.0 }), { 1, 1, 1 }
        );
        BOOST_CHECK_THROW((void)GridIndexer(*non_finite), queso::Exception);
    }

    BOOST_AUTO_TEST_CASE(RejectsOneThinCellDirection)
    {
        auto p_settings = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 1e-15, 1.0, 1.0 }), { 1, 1, 1 });
        BOOST_CHECK_THROW((void)GridIndexer(*p_settings), queso::Exception);
    }

    BOOST_AUTO_TEST_CASE(NonBinaryGridPlanesPreserveOwnership)
    {
        auto ownership_settings = MakeGridSettings(MakeBox({ 0.1, 0.1, 0.1 }, { 1.1, 1.1, 1.1 }), { 4, 4, 4 });
        const GridIndexer ownership_indexer(*ownership_settings);
        const PointType internal_planes{
            ownership_indexer.GetFaceCoordinate(ownership_indexer.GetFace(0, GridIndexer::Direction::x_forward)),
            ownership_indexer.GetFaceCoordinate(ownership_indexer.GetFace(0, GridIndexer::Direction::y_forward)),
            ownership_indexer.GetFaceCoordinate(ownership_indexer.GetFace(0, GridIndexer::Direction::z_forward))
        };
        QuESo_CHECK_EQUAL(
            ownership_indexer.GetMatrixIndicesFromVectorIndex(ownership_indexer.GetContainingCell(internal_planes)),
            (Vector3i{ 1, 1, 1 })
        );
        QuESo_CHECK_EQUAL(ownership_indexer.GetContainingCell(PointType{ 0.1, 0.1, 0.1 }), IndexType{ 0 });
        QuESo_CHECK_EQUAL(
            ownership_indexer.GetMatrixIndicesFromVectorIndex(
                ownership_indexer.GetContainingCell(PointType{ 1.1, 1.1, 1.1 })
            ),
            (Vector3i{ 3, 3, 3 })
        );
    }

    BOOST_AUTO_TEST_CASE(PhysicalFaceIdsAreDenseOnAsymmetricGrid)
    {
        constexpr Vector3i counts{ 2, 3, 4 };
        auto p_settings = MakeGridSettings(MakeBox({ -1.0, 0.5, 2.0 }, { 3.0, 6.5, 10.0 }), counts);
        const GridIndexer grid_indexer(*p_settings);
        std::set<IndexType> face_ids;
        for (IndexType cell = 0; cell < grid_indexer.NumberOfElements(); ++cell) {
            for (const auto direction : EnumRange<GridIndexer::Direction>()) {
                const GridFaceId face = grid_indexer.GetFace(cell, direction);
                QuESo_CHECK(face.value < grid_indexer.NumberOfFaces());
                face_ids.insert(face.value);
                QuESo_CHECK(grid_indexer.GetAdjacentCell(face, grid_indexer.GetAdjacentSide(face, cell)) == cell);
                QuESo_CHECK_EQUAL(grid_indexer.GetFaceAxis(face), static_cast<IndexType>(direction) / 2);
                if (!grid_indexer.IsEnd(cell, direction)) {
                    const auto [neighbor, info] = grid_indexer.GetNextIndex(cell, direction);
                    QuESo_CHECK(info == GridIndexer::IndexInfo::middle);
                    QuESo_CHECK(grid_indexer.GetFace(neighbor, GridIndexer::ReverseDirection(direction)) == face);
                }
            }
        }
        QuESo_CHECK_EQUAL(face_ids.size(), grid_indexer.NumberOfFaces());
    }

    BOOST_AUTO_TEST_CASE(PhysicalFaceNavigation)
    {
        auto p_settings = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 2.0, 1.0, 1.0 }), { 2, 1, 1 });
        const GridIndexer grid_indexer(*p_settings);

        QuESo_CHECK_EQUAL(grid_indexer.NumberOfFaces(), 11);
        const GridFaceId shared_from_left = grid_indexer.GetFace(0, GridIndexer::Direction::x_forward);
        const GridFaceId shared_from_right = grid_indexer.GetFace(1, GridIndexer::Direction::x_backward);
        QuESo_CHECK(shared_from_left == shared_from_right);
        QuESo_CHECK_EQUAL(*grid_indexer.GetAdjacentCell(shared_from_left, GridFaceSide::negative), IndexType{ 0 });
        QuESo_CHECK_EQUAL(*grid_indexer.GetAdjacentCell(shared_from_left, GridFaceSide::positive), IndexType{ 1 });
        QuESo_CHECK_EQUAL(grid_indexer.GetCanonicalAdjacentCell(shared_from_left), IndexType{ 1 });

        const GridFaceId negative_exterior = grid_indexer.GetFace(0, GridIndexer::Direction::x_backward);
        QuESo_CHECK(!grid_indexer.GetAdjacentCell(negative_exterior, GridFaceSide::negative));
        QuESo_CHECK_EQUAL(*grid_indexer.GetAdjacentCell(negative_exterior, GridFaceSide::positive), IndexType{ 0 });
        QuESo_CHECK_EQUAL(grid_indexer.GetCanonicalAdjacentCell(negative_exterior), IndexType{ 0 });

        const GridFaceId positive_exterior = grid_indexer.GetFace(1, GridIndexer::Direction::x_forward);
        QuESo_CHECK_EQUAL(*grid_indexer.GetAdjacentCell(positive_exterior, GridFaceSide::negative), IndexType{ 1 });
        QuESo_CHECK(!grid_indexer.GetAdjacentCell(positive_exterior, GridFaceSide::positive));
        QuESo_CHECK_EQUAL(grid_indexer.GetCanonicalAdjacentCell(positive_exterior), IndexType{ 1 });
        QuESo_CHECK_EQUAL(Opposite(GridFaceSide::negative), GridFaceSide::positive);
    }

    BOOST_AUTO_TEST_SUITE_END()

}  // End namespace Testing
}  // End namespace queso
