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

    Settings settings;
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::grid_type, BackgroundGridType::b_spline_grid);
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.first);
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.second);
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.first);
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.second);
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);

    GridIndexer grid_indexer(settings);
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

    Settings settings;
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::grid_type, BackgroundGridType::hexahedral_fe_grid);
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, bounds_xyz.first);
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, bounds_xyz.second);
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, bounds_uvw.first);
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, bounds_uvw.second);
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);

    GridIndexer grid_indexer(settings);
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

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso