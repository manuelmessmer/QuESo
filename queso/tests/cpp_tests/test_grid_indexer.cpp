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
#include "queso/containers/grid_indexer.hpp"
#include "queso/utilities/math_utilities.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( GridIndexerTestSuite )

BOOST_AUTO_TEST_CASE(GridIndexerBSplineMeshTest) {
    QuESo_INFO << "Testing :: Test Grid Indexer :: BSpline Mesh" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-1.0, -0.5, 1.0}, {5.0, 10.5, 13.0} );
    const BoundingBoxType bounds_uvw = MakeBox( {-1.0, -0.5, 1.0}, {5.0, 10.5, 13.0} );

    const Vector3i number_of_elements{5, 10, 7};

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.second);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.second);
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);

    GridIndexer grid_indexer(r_settings);
    const auto delta = Math::Subtract( bounds_xyz.second,  bounds_xyz.first );
    double volume = 0.0;
    for( IndexType i = 0; i < number_of_elements[0]; ++i){
        for( IndexType j = 0; j < number_of_elements[1]; ++j){
            for( IndexType k = 0; k < number_of_elements[2]; ++k){
                IndexType index = grid_indexer.GetVectorIndexFromMatrixIndices(i, j, k);
                auto indices = grid_indexer.GetMatrixIndicesFromVectorIndex(index);
                QuESo_CHECK_EQUAL(i, indices[0]);
                QuESo_CHECK_EQUAL(j, indices[1]);
                QuESo_CHECK_EQUAL(k, indices[2]);

                auto box_1_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(index);
                auto box_2_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(i, j, k);
                auto box_3_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(indices);

                QuESo_CHECK_POINT_NEAR( box_1_xyz.first, box_2_xyz.first, 1e-12 );
                QuESo_CHECK_POINT_NEAR( box_1_xyz.first, box_3_xyz.first, 1e-12 );
                QuESo_CHECK_POINT_NEAR( box_1_xyz.second, box_2_xyz.second, 1e-12 );
                QuESo_CHECK_POINT_NEAR( box_1_xyz.second, box_3_xyz.second, 1e-12 );

                auto box_1_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(index);
                auto box_2_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(i, j, k);
                auto box_3_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(indices);

                QuESo_CHECK_POINT_NEAR( box_1_uvw.first, box_1_xyz.first, 1e-12 );
                QuESo_CHECK_POINT_NEAR( box_1_uvw.second, box_1_xyz.second, 1e-12 );
                QuESo_CHECK_POINT_NEAR( box_2_uvw.first, box_2_xyz.first, 1e-12 );
                QuESo_CHECK_POINT_NEAR( box_2_uvw.second, box_2_xyz.second, 1e-12 );
                QuESo_CHECK_POINT_NEAR( box_3_uvw.first, box_3_xyz.first, 1e-12 );
                QuESo_CHECK_POINT_NEAR( box_3_uvw.second, box_3_xyz.second, 1e-12 );

                auto delta_box = Math::Subtract(box_1_xyz.second, box_1_xyz.first);
                QuESo_CHECK_LT( std::abs(delta[0]/number_of_elements[0] - delta_box[0]) / delta_box[0], 1e-12);
                QuESo_CHECK_LT( std::abs(delta[1]/number_of_elements[1] - delta_box[1]) / delta_box[1], 1e-12);
                QuESo_CHECK_LT( std::abs(delta[2]/number_of_elements[2] - delta_box[2]) / delta_box[2], 1e-12);

                volume += delta_box[0]*delta_box[1]*delta_box[2];

                IndexType index_1 = grid_indexer.GetVectorIndexFromMatrixIndices(i, j, k);
                IndexType index_2 = grid_indexer.GetVectorIndexFromMatrixIndices(indices);
                QuESo_CHECK_EQUAL(index_1, index);
                QuESo_CHECK_EQUAL(index_2, index);
            }
        }
    }

    double volume_ref = delta[0]*delta[1]*delta[2];
    QuESo_CHECK_LT( std::abs(volume - volume_ref) / volume_ref, 1e-12);
}


BOOST_AUTO_TEST_CASE(GridIndexerFEMeshTest) {
    QuESo_INFO << "Testing :: Test Grid Indexer :: Hexahedral FE Mesh" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-1.0, -0.5, 1.0}, {5.0, 10.5, 13.0} );
    const BoundingBoxType bounds_uvw = MakeBox( {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0} );

    const Vector3i number_of_elements{5, 10, 7};

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.second);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.second);
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);

    GridIndexer grid_indexer(r_settings);
    const auto delta = Math::Subtract( bounds_xyz.second,  bounds_xyz.first );
    double volume = 0.0;
    for( IndexType i = 0; i < number_of_elements[0]; ++i){
        for( IndexType j = 0; j < number_of_elements[1]; ++j){
            for( IndexType k = 0; k < number_of_elements[2]; ++k){
                IndexType index = grid_indexer.GetVectorIndexFromMatrixIndices(i, j, k);
                auto indices = grid_indexer.GetMatrixIndicesFromVectorIndex(index);
                QuESo_CHECK_EQUAL(i, indices[0]);
                QuESo_CHECK_EQUAL(j, indices[1]);
                QuESo_CHECK_EQUAL(k, indices[2]);

                auto box_1_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(index);
                auto box_2_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(i, j, k);
                auto box_3_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(indices);

                QuESo_CHECK_POINT_NEAR(box_1_xyz.first, box_2_xyz.first, 1e-12 );
                QuESo_CHECK_POINT_NEAR(box_1_xyz.first, box_3_xyz.first, 1e-12 );
                QuESo_CHECK_POINT_NEAR(box_1_xyz.second, box_2_xyz.second, 1e-12 );
                QuESo_CHECK_POINT_NEAR(box_1_xyz.second, box_3_xyz.second, 1e-12 );

                auto box_1_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(index);
                auto box_2_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(i, j, k);
                auto box_3_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(indices);

                QuESo_CHECK_POINT_NEAR(box_1_uvw.first, PointType({-1.0, -1.0, -1.0}), 1e-12);
                QuESo_CHECK_POINT_NEAR(box_1_uvw.second, PointType({1.0, 1.0, 1.0}), 1e-12);
                QuESo_CHECK_POINT_NEAR(box_2_uvw.first, PointType({-1.0, -1.0, -1.0}), 1e-12);
                QuESo_CHECK_POINT_NEAR(box_2_uvw.second, PointType({1.0, 1.0, 1.0}), 1e-12);
                QuESo_CHECK_POINT_NEAR(box_3_uvw.first, PointType({-1.0, -1.0, -1.0}), 1e-12);
                QuESo_CHECK_POINT_NEAR(box_3_uvw.second, PointType({1.0, 1.0, 1.0}), 1e-12);

                auto delta_box = Math::Subtract(box_1_xyz.second, box_1_xyz.first);
                QuESo_CHECK_LT( std::abs(delta[0]/number_of_elements[0] - delta_box[0]) / delta_box[0], 1e-12);
                QuESo_CHECK_LT( std::abs(delta[1]/number_of_elements[1] - delta_box[1]) / delta_box[1], 1e-12);
                QuESo_CHECK_LT( std::abs(delta[2]/number_of_elements[2] - delta_box[2]) / delta_box[2], 1e-12);

                volume += delta_box[0]*delta_box[1]*delta_box[2];

                IndexType index_1 = grid_indexer.GetVectorIndexFromMatrixIndices(i, j, k);
                IndexType index_2 = grid_indexer.GetVectorIndexFromMatrixIndices(indices);
                QuESo_CHECK_EQUAL(index_1, index);
                QuESo_CHECK_EQUAL(index_2, index);
            }
        }
    }

    double volume_ref = delta[0]*delta[1]*delta[2];
    QuESo_CHECK_LT( std::abs(volume - volume_ref) / volume_ref, 1e-12);
}

bool contains(std::vector<IndexType>& v, IndexType test_value){
    if(std::find(v.begin(), v.end(), test_value) != v.end()) {
        return true;
    }
    return false;
}

BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingGlobalXTest) {
    QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Global Partition X" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-1.0, -0.5, 1.0}, {5.0, 10.5, 13.0} );
    const BoundingBoxType bounds_uvw = MakeBox( {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0} );

    const Vector3i number_of_elements{3, 4, 5};

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.second);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.second);
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);

    GridIndexer grid_indexer(r_settings);

    // Walk foward
    IndexType index = 0;
    GridIndexer::IndexInfo index_info{};
    IndexType i = 0;
    while(index < grid_indexer.NumberOfElements()-1 ){
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if ((index+1)%3 == 0) { index_info_ref = GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetNextIndexX(index);
        ++i;
        QuESo_CHECK_EQUAL(index, i);
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
    }
    QuESo_CHECK_EQUAL(index, 59);
    std::tie(index, index_info) = grid_indexer.GetNextIndexX(index);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end );
    QuESo_CHECK_EQUAL(index, 59);
    --i;

    // Walk backwards
    while( index > 0) {
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if ((index)%3 == 0) { index_info_ref = GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetPreviousIndexX(index);
        QuESo_CHECK_EQUAL(index, i);
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
        --i;
    }
    QuESo_CHECK_EQUAL(index, 0);
    std::tie(index, index_info) = grid_indexer.GetPreviousIndexX(index);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
    QuESo_CHECK_EQUAL(index, 0);
}

BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingGlobalYTest) {
    QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Global Partition Y" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-1.0, -0.5, 1.0}, {5.0, 10.5, 13.0} );
    const BoundingBoxType bounds_uvw = MakeBox( {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0} );

    const Vector3i number_of_elements{3, 4, 5};

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.second);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.second);
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);

    GridIndexer grid_indexer(r_settings);
    // Walk foward
    IndexType index = 0;
    GridIndexer::IndexInfo index_info;
    IndexType i = 1;
    std::vector<IndexType> order_y = {0,3,6,9,1,4,7,10,2,5,8,11,12,15,18,21,13,16,19,22,14,17,20,23,24,27,30,
        33,25,28,31,34,26,29,32,35,36,39,42,45,37,40,43,46,38,41,44,47,48,51,54,57,49,52,55,58,50,53,56,59};
    std::vector<IndexType> local_ends_y = {9, 10, 11, 21, 22, 23, 33, 34, 35, 45, 46, 47, 57, 58};
    while(index < grid_indexer.NumberOfElements()-1 ){
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends_y, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetNextIndexY(index);
        QuESo_CHECK_EQUAL(order_y[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
        ++i;
    }

    QuESo_CHECK_EQUAL(index, 59);
    std::tie(index, index_info) = grid_indexer.GetNextIndexY(index);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end );
    QuESo_CHECK_EQUAL(index, 59);
    --i;

    local_ends_y = {0, 1, 2, 12, 13, 14, 24, 25, 26, 36, 37, 38, 48, 49, 50};
    while( index > 0) {
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends_y, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetPreviousIndexY(index);
        --i;
        QuESo_CHECK_EQUAL(order_y[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
    }
    QuESo_CHECK_EQUAL(index, 0);
    std::tie(index, index_info) = grid_indexer.GetPreviousIndexZ(index);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
    QuESo_CHECK_EQUAL(index, 0);
}


BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingGlobalZTest) {
    QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Global Partition Z" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-1.0, -0.5, 1.0}, {5.0, 10.5, 13.0} );
    const BoundingBoxType bounds_uvw = MakeBox( {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0} );

    const Vector3i number_of_elements{3, 4, 5};

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.second);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.second);
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);

    GridIndexer grid_indexer(r_settings);
    // Walk foward
    IndexType index = 0;
    GridIndexer::IndexInfo index_info;
    IndexType i = 1;
    std::vector<IndexType> order_z = {0,12,24,36,48,1,13,25,37,49,2,14,26,38,50,3,15,27,39,51,4,16,28,40,52,5,17,29,41,53,
        6,18,30,42,54,7,19,31,43,55,8,20,32,44,56,9,21,33,45,57,10,22,34,46,58,11,23,35,47,59};

    std::vector<IndexType> local_ends_z = {48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58};
    while(index < grid_indexer.NumberOfElements()-1  ){
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends_z, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetNextIndexZ(index);
        QuESo_CHECK_EQUAL(order_z[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
        ++i;
    }

    QuESo_CHECK_EQUAL(index, 59);
    std::tie(index, index_info) = grid_indexer.GetNextIndexZ(index);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end );
    QuESo_CHECK_EQUAL(index, 59);
    --i;

    // Walk backwards
    local_ends_z = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    while( index > 0) {
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends_z, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetPreviousIndexZ(index);
        --i;
        QuESo_CHECK_EQUAL(order_z[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
    }
    QuESo_CHECK_EQUAL(index, 0);
    std::tie(index, index_info) = grid_indexer.GetPreviousIndexZ(index);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
    QuESo_CHECK_EQUAL(index, 0);
}

BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingLocalXTest) {
    QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Local Partition X" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-1.0, -0.5, 1.0}, {5.0, 10.5, 13.0} );
    const BoundingBoxType bounds_uvw = MakeBox( {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0} );

    const Vector3i number_of_elements{3, 4, 5};

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.second);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.second);
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);

    GridIndexer grid_indexer(r_settings);

    PartitionBoxType local_partition_1 = std::make_pair(Vector3i({1, 1, 1}), Vector3i({2, 3, 3}));
    // Walk foward
    IndexType index = 16;
    GridIndexer::IndexInfo index_info;
    IndexType i = 1;
    std::vector<IndexType> order = {16, 17, 19, 20, 22, 23, 28, 29, 31, 32, 34, 35, 40, 41, 43, 44, 46, 47};
    std::vector<IndexType> local_ends = {17, 20, 23, 29, 32, 35, 41, 44, 47};
    while(index < 47 ){
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetNextIndexX(index, local_partition_1);
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
        ++i;
    }

    QuESo_CHECK_EQUAL(index, 47);
    std::tie(index, index_info) = grid_indexer.GetNextIndexX(index, local_partition_1);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end );
    QuESo_CHECK_EQUAL(index, 47);
    --i;

    // Walk backwards
    local_ends = {16, 19, 22, 28, 31, 34, 40, 43, 46};
    while( index > 16) {
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetPreviousIndexX(index, local_partition_1);
        --i;
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
    }
    QuESo_CHECK_EQUAL(index, 16);
    std::tie(index, index_info) = grid_indexer.GetPreviousIndexX(index, local_partition_1);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
    QuESo_CHECK_EQUAL(index, 16);

    PartitionBoxType local_partition_2 = std::make_pair(Vector3i({0, 0, 0}), Vector3i({1, 1, 2}));
    // Walk foward
    index = 0;
    i = 1;
    order = {0, 1, 3, 4, 12, 13, 15, 16, 24, 25, 27, 28};
    local_ends = {1, 4, 13, 16, 25, 28};
    while(index < 28 ){
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetNextIndexX(index, local_partition_2);
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
        ++i;
    }

    QuESo_CHECK_EQUAL(index, 28);
    std::tie(index, index_info) = grid_indexer.GetNextIndexX(index, local_partition_2);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end );
    QuESo_CHECK_EQUAL(index, 28);
    --i;

    // Walk backwards
    local_ends = {0, 3, 12, 15, 24, 27};
    while( index > 0) {
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetPreviousIndexX(index, local_partition_2);
        --i;
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
    }
    QuESo_CHECK_EQUAL(index, 0);
    std::tie(index, index_info) = grid_indexer.GetPreviousIndexX(index, local_partition_2);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
    QuESo_CHECK_EQUAL(index, 0);
}

BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingLocalYTest) {
    QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Local Partition Y" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-1.0, -0.5, 1.0}, {5.0, 10.5, 13.0} );
    const BoundingBoxType bounds_uvw = MakeBox( {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0} );

    const Vector3i number_of_elements{3, 4, 5};

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.second);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.second);
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);

    GridIndexer grid_indexer(r_settings);

    PartitionBoxType local_partition_1 = std::make_pair(Vector3i({1, 1, 1}), Vector3i({2, 3, 3}));
    // Walk foward
    IndexType index = 16;
    GridIndexer::IndexInfo index_info;
    IndexType i = 1;
    std::vector<IndexType> order = {16, 19, 22, 17, 20, 23, 28, 31, 34, 29, 32, 35, 40, 43, 46, 41, 44, 47};
    std::vector<IndexType> local_ends = {22, 23, 34, 35, 46, 47};
    while(index < 47 ){
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetNextIndexY(index, local_partition_1);
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
        ++i;
    }

    QuESo_CHECK_EQUAL(index, 47);
    std::tie(index, index_info) = grid_indexer.GetNextIndexY(index, local_partition_1);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end );
    QuESo_CHECK_EQUAL(index, 47);
    --i;

    // Walk backwards
    local_ends = {16, 17, 28, 29, 40, 41};
    while( index > 16) {
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetPreviousIndexY(index, local_partition_1);
        --i;
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
    }
    QuESo_CHECK_EQUAL(index, 16);
    std::tie(index, index_info) = grid_indexer.GetPreviousIndexY(index, local_partition_1);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
    QuESo_CHECK_EQUAL(index, 16);


    PartitionBoxType local_partition_2 = std::make_pair(Vector3i({0, 0, 0}), Vector3i({1, 1, 2}));
    // Walk foward
    index = 0;
    i = 1;
    order = {0, 3, 1, 4, 12, 15, 13, 16, 24, 27, 25, 28};
    local_ends = {3, 4, 15, 16, 27, 28};
    while(index < 28 ){
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetNextIndexY(index, local_partition_2);
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
        ++i;
    }

    QuESo_CHECK_EQUAL(index, 28);
    std::tie(index, index_info) = grid_indexer.GetNextIndexY(index, local_partition_2);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end );
    QuESo_CHECK_EQUAL(index, 28);
    --i;

    // Walk backwards
    local_ends = {0, 1, 12, 13, 24, 25};
    while( index > 0) {
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetPreviousIndexY(index, local_partition_2);
        --i;
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
    }
    QuESo_CHECK_EQUAL(index, 0);
    std::tie(index, index_info) = grid_indexer.GetPreviousIndexY(index, local_partition_2);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
    QuESo_CHECK_EQUAL(index, 0);
}


BOOST_AUTO_TEST_CASE(GridIndexerIndexWalkingLocalZTest) {
    QuESo_INFO << "Testing :: Test Grid Indexer :: Test Walk Through Local Partition Z" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-1.0, -0.5, 1.0}, {5.0, 10.5, 13.0} );
    const BoundingBoxType bounds_uvw = MakeBox( {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0} );

    const Vector3i number_of_elements{3, 4, 5};

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.second);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.first);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.second);
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);

    GridIndexer grid_indexer(r_settings);

    PartitionBoxType local_partition_1 = std::make_pair(Vector3i({1, 1, 1}), Vector3i({2, 3, 3}));
    // Walk foward
    IndexType index = 16;
    GridIndexer::IndexInfo index_info;
    IndexType i = 1;
    std::vector<IndexType> order = {16, 28, 40, 17, 29, 41, 19, 31, 43, 20, 32, 44, 22, 34, 46, 23, 35, 47};
    std::vector<IndexType> local_ends = {40, 41, 43, 44, 46, 47};
    while(index < 47 ){
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetNextIndexZ(index, local_partition_1);
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
        ++i;
    }

    QuESo_CHECK_EQUAL(index, 47);
    std::tie(index, index_info) = grid_indexer.GetNextIndexZ(index, local_partition_1);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end );
    QuESo_CHECK_EQUAL(index, 47);
    --i;

    // Walk backwards
    local_ends = {16, 17, 19, 20, 22, 23};
    while( index > 16) {
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetPreviousIndexZ(index, local_partition_1);
        --i;
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
    }
    QuESo_CHECK_EQUAL(index, 16);
    std::tie(index, index_info) = grid_indexer.GetPreviousIndexZ(index, local_partition_1);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
    QuESo_CHECK_EQUAL(index, 16);

    PartitionBoxType local_partition_2 = std::make_pair(Vector3i({0, 0, 0}), Vector3i({1, 1, 2}));
    // Walk foward
    index = 0;
    i = 1;
    order = {0, 12, 24, 1, 13, 25, 3, 15, 27, 4, 16, 28};
    local_ends = {27, 28, 24, 25};
    while(index < 28 ){
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetNextIndexZ(index, local_partition_2);
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
        ++i;
    }

    QuESo_CHECK_EQUAL(index, 28);
    std::tie(index, index_info) = grid_indexer.GetNextIndexZ(index, local_partition_2);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end );
    QuESo_CHECK_EQUAL(index, 28);
    --i;

    // Walk backwards
    local_ends = {0, 1, 3, 4};
    while( index > 0) {
        GridIndexer::IndexInfo index_info_ref = GridIndexer::IndexInfo::middle;
        if( contains(local_ends, index) ) { index_info_ref =  GridIndexer::IndexInfo::local_end; }
        std::tie(index, index_info) = grid_indexer.GetPreviousIndexZ(index, local_partition_2);
        --i;
        QuESo_CHECK_EQUAL(order[i], index)
        QuESo_CHECK_EQUAL(index_info, index_info_ref);
    }
    QuESo_CHECK_EQUAL(index, 0);
    std::tie(index, index_info) = grid_indexer.GetPreviousIndexZ(index, local_partition_2);
    QuESo_CHECK_EQUAL(index_info, GridIndexer::IndexInfo::global_end);
    QuESo_CHECK_EQUAL(index, 0);
}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso