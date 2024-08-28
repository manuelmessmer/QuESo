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
#include "queso/utilities/mapping_utilities.h"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( VoxelIndexerTestSuite )

BOOST_AUTO_TEST_CASE(VoxelIndexerTest) {
    QuESo_INFO << "Testing :: Test VoxelIndexing Utilities :: VoxelIndexer" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-1.0, -0.5, 1.0}, {5.0, 10.5, 13.0} );
    const BoundingBoxType bounds_uvw = MakeBox( { -10.0,  5.0, 2.0}, {-2.0, 7.0, 100.0} );

    const Vector3i number_of_elements{5, 10, 7};

    Parameters parameters( {Component("lower_bound_xyz", bounds_xyz.first),
                            Component("upper_bound_xyz", bounds_xyz.second),
                            Component("lower_bound_uvw", bounds_uvw.first),
                            Component("upper_bound_uvw", bounds_uvw.second),
                            Component("number_of_elements", number_of_elements) } );

    VoxelIndexer voxel_indexer(parameters);
    const auto delta = Math::Subtract( bounds_xyz.second,  bounds_xyz.first );
    double volume = 0.0;
    for( IndexType i = 0; i < number_of_elements[0]; ++i){
        for( IndexType j = 0; j < number_of_elements[1]; ++j){
            for( IndexType k = 0; k < number_of_elements[2]; ++k){
                IndexType index = voxel_indexer.GetVectorIndexFromMatrixIndices(i, j, k);
                auto indices = voxel_indexer.GetMatrixIndicesFromVectorIndex(index);
                QuESo_CHECK_EQUAL(i, indices[0]);
                QuESo_CHECK_EQUAL(j, indices[1]);
                QuESo_CHECK_EQUAL(k, indices[2]);

                auto box_1 = voxel_indexer.GetBoundingBoxXYZFromIndex(index);
                auto box_2 = voxel_indexer.GetBoundingBoxXYZFromIndex(i, j, k);
                auto box_3 = voxel_indexer.GetBoundingBoxXYZFromIndex(indices);

                QuESo_CHECK_LT( Math::Norm( Math::Subtract(box_1.first,  box_2.first) ), 1e-12 );
                QuESo_CHECK_LT( Math::Norm( Math::Subtract(box_1.first, box_3.first) ), 1e-12 );
                QuESo_CHECK_LT( Math::Norm( Math::Subtract(box_1.second, box_2.second) ), 1e-12 );
                QuESo_CHECK_LT( Math::Norm( Math::Subtract(box_1.second, box_3.second) ), 1e-12 );

                auto delta_box = Math::Subtract(box_1.second, box_1.first);
                QuESo_CHECK_LT( std::abs(delta[0]/number_of_elements[0] - delta_box[0]) / delta_box[0], 1e-12);
                QuESo_CHECK_LT( std::abs(delta[1]/number_of_elements[1] - delta_box[1]) / delta_box[1], 1e-12);
                QuESo_CHECK_LT( std::abs(delta[2]/number_of_elements[2] - delta_box[2]) / delta_box[2], 1e-12);

                volume += delta_box[0]*delta_box[1]*delta_box[2];

                IndexType index_1 = voxel_indexer.GetVectorIndexFromMatrixIndices(i, j, k);
                IndexType index_2 = voxel_indexer.GetVectorIndexFromMatrixIndices(indices);
                QuESo_CHECK_EQUAL(index_1, index);
                QuESo_CHECK_EQUAL(index_2, index);
            }
        }
    }

    double volume_ref = delta[0]*delta[1]*delta[2];
    QuESo_CHECK_LT( std::abs(volume - volume_ref) / volume_ref, 1e-12);
}

BOOST_AUTO_TEST_CASE(VoxelIndexingTest) {
    QuESo_INFO << "Testing :: Test VoxelIndexing Utilities :: VoxelIndexing" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-25, -111.44, 7.89}, {78.67, -35.68, 18.99} );
    const BoundingBoxType bounds_uvw = MakeBox( { -10.0,  -2.2, 2.0}, {2.0, 10.0, 17.0} );

    const Vector3i number_of_elements{15, 11, 17};

    Parameters parameters( {Component("lower_bound_xyz", bounds_xyz.first),
                            Component("upper_bound_xyz", bounds_xyz.second),
                            Component("lower_bound_uvw", bounds_uvw.first),
                            Component("upper_bound_uvw", bounds_uvw.second),
                            Component("number_of_elements", number_of_elements) } );

    const auto delta = Math::Subtract( bounds_xyz.second, bounds_xyz.first );
    double volume = 0.0;
    for( IndexType i = 0; i < number_of_elements[0]; ++i){
        for( IndexType j = 0; j < number_of_elements[1]; ++j){
            for( IndexType k = 0; k < number_of_elements[2]; ++k){
                IndexType index = VoxelIndexing::GetVectorIndexFromMatrixIndices(i, j, k, number_of_elements);
                auto indices = VoxelIndexing::GetMatrixIndicesFromVectorIndex(index, number_of_elements);
                QuESo_CHECK_EQUAL(i, indices[0]);
                QuESo_CHECK_EQUAL(j, indices[1]);
                QuESo_CHECK_EQUAL(k, indices[2]);

                auto box_1 = VoxelIndexing::GetBoundingBoxFromIndex(index, bounds_xyz.first, bounds_xyz.second, number_of_elements);
                auto box_2 = VoxelIndexing::GetBoundingBoxFromIndex(i, j, k, bounds_xyz.first, bounds_xyz.second, number_of_elements);
                auto box_3 = VoxelIndexing::GetBoundingBoxFromIndex(indices, bounds_xyz.first, bounds_xyz.second, number_of_elements);

                QuESo_CHECK_LT( Math::Norm( Math::Subtract(box_1.first, box_2.first) ), 1e-12 );
                QuESo_CHECK_LT( Math::Norm( Math::Subtract(box_1.first, box_3.first) ), 1e-12 );
                QuESo_CHECK_LT( Math::Norm( Math::Subtract(box_1.second, box_2.second) ), 1e-12 );
                QuESo_CHECK_LT( Math::Norm( Math::Subtract(box_1.second, box_3.second) ), 1e-12 );

                auto delta_box = Math::Subtract(box_1.second, box_1.first);
                QuESo_CHECK_LT( std::abs(delta[0]/number_of_elements[0] - delta_box[0]) / delta_box[0], 1e-12);
                QuESo_CHECK_LT( std::abs(delta[1]/number_of_elements[1] - delta_box[1]) / delta_box[1], 1e-12);
                QuESo_CHECK_LT( std::abs(delta[2]/number_of_elements[2] - delta_box[2]) / delta_box[2], 1e-12);

                volume += delta_box[0]*delta_box[1]*delta_box[2];

                IndexType index_1 = VoxelIndexing::GetVectorIndexFromMatrixIndices(i, j, k, number_of_elements);
                IndexType index_2 = VoxelIndexing::GetVectorIndexFromMatrixIndices(indices, number_of_elements);
                QuESo_CHECK_EQUAL(index_1, index);
                QuESo_CHECK_EQUAL(index_2, index);

                PointType point_global{i+100.0, j-44.27, k+88.90};
                PointType point_param =  VoxelIndexing::PointFromGlobalToParam(point_global, bounds_xyz, bounds_uvw);
                PointType point_global_2 =  VoxelIndexing::PointFromParamToGlobal(point_param, bounds_xyz, bounds_uvw);

                QuESo_CHECK_LT( Math::Norm( Math::Subtract(point_global, point_global_2) ), 1e-12 );
            }
        }
    }

    double volume_ref = delta[0]*delta[1]*delta[2];
    QuESo_CHECK_LT( std::abs(volume - volume_ref) / volume_ref , 1e-12);
}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso