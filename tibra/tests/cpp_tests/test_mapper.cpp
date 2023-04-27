// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "utilities/mapping_utilities.h"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( MapperTestSuite )

BOOST_AUTO_TEST_CASE(MapperTest) {
    TIBRA_INFO << "Testing :: Test Mapping Utilities :: Mapper" << std::endl;

    const PointType lower_bound(-1.0, -0.5, 1.0);
    const PointType upper_bound(5.0, 10.5, 13.0);
    const Vector3i number_of_elements(5, 10, 7);

    Parameters parameters( {Component("lower_bound", lower_bound),
                        Component("upper_bound", upper_bound),
                        Component("number_of_elements", number_of_elements) } );

    Mapper mapper(parameters);
    const auto delta = upper_bound - lower_bound;
    double volume = 0.0;
    for( IndexType i = 0; i < number_of_elements[0]; ++i){
        for( IndexType j = 0; j < number_of_elements[1]; ++j){
            for( IndexType k = 0; k < number_of_elements[2]; ++k){
                IndexType index = mapper.GetVectorIndexFromMatrixIndices(i, j, k);
                auto indices = mapper.GetMatrixIndicesFromVectorIndex(index);
                BOOST_CHECK_EQUAL(i, indices[0]);
                BOOST_CHECK_EQUAL(j, indices[1]);
                BOOST_CHECK_EQUAL(k, indices[2]);

                auto box_1 = mapper.GetBoundingBoxFromIndex(index);
                auto box_2 = mapper.GetBoundingBoxFromIndex(i, j, k);
                auto box_3 = mapper.GetBoundingBoxFromIndex(indices);

                BOOST_CHECK_LT( (box_1.first - box_2.first).Norm(), 1e-12 );
                BOOST_CHECK_LT( (box_1.first - box_3.first).Norm(), 1e-12 );
                BOOST_CHECK_LT( (box_1.second - box_2.second).Norm(), 1e-12 );
                BOOST_CHECK_LT( (box_1.second - box_3.second).Norm(), 1e-12 );

                auto delta_box = box_1.second - box_1.first;
                BOOST_CHECK_LT( std::abs(delta[0]/number_of_elements[0] - delta_box[0]) / delta_box[0], 1e-12);
                BOOST_CHECK_LT( std::abs(delta[1]/number_of_elements[1] - delta_box[1]) / delta_box[1], 1e-12);
                BOOST_CHECK_LT( std::abs(delta[2]/number_of_elements[2] - delta_box[2]) / delta_box[2], 1e-12);

                volume += delta_box[0]*delta_box[1]*delta_box[2];

                IndexType index_1 = mapper.GetVectorIndexFromMatrixIndices(i, j, k);
                IndexType index_2 = mapper.GetVectorIndexFromMatrixIndices(indices);
                BOOST_CHECK_EQUAL(index_1, index);
                BOOST_CHECK_EQUAL(index_2, index);

                PointType point_global(i+0.1232321, j-123.44, k+0.333);
                PointType point_param =  mapper.PointFromGlobalToParam(point_global);
                PointType point_global_2 =  mapper.PointFromParamToGlobal(point_param);

                BOOST_CHECK_LT( (point_global - point_global_2).Norm(), 1e-12 );
            }
        }
    }

    double volume_ref = delta[0]*delta[1]*delta[2];
    BOOST_CHECK_LT( std::abs(volume - volume_ref) / volume_ref, 1e-12);
}

BOOST_AUTO_TEST_CASE(MappingTest) {
    TIBRA_INFO << "Testing :: Test Mapping Utilities :: Mapping" << std::endl;

    const PointType lower_bound(-25, -111.44, 7.89);
    const PointType upper_bound(78.67, -35.68, 18.99);
    const Vector3i number_of_elements(15, 11, 17);

    Parameters parameters( {Component("lower_bound", lower_bound),
                        Component("upper_bound", upper_bound),
                        Component("number_of_elements", number_of_elements) } );

    const auto delta = upper_bound - lower_bound;
    double volume = 0.0;
    for( IndexType i = 0; i < number_of_elements[0]; ++i){
        for( IndexType j = 0; j < number_of_elements[1]; ++j){
            for( IndexType k = 0; k < number_of_elements[2]; ++k){
                IndexType index = Mapping::GetVectorIndexFromMatrixIndices(i, j, k, number_of_elements);
                auto indices = Mapping::GetMatrixIndicesFromVectorIndex(index, number_of_elements);
                BOOST_CHECK_EQUAL(i, indices[0]);
                BOOST_CHECK_EQUAL(j, indices[1]);
                BOOST_CHECK_EQUAL(k, indices[2]);

                auto box_1 = Mapping::GetBoundingBoxFromIndex(index, lower_bound, upper_bound, number_of_elements);
                auto box_2 = Mapping::GetBoundingBoxFromIndex(i, j, k, lower_bound, upper_bound, number_of_elements);
                auto box_3 = Mapping::GetBoundingBoxFromIndex(indices, lower_bound, upper_bound, number_of_elements);

                BOOST_CHECK_LT( (box_1.first - box_2.first).Norm(), 1e-12 );
                BOOST_CHECK_LT( (box_1.first - box_3.first).Norm(), 1e-12 );
                BOOST_CHECK_LT( (box_1.second - box_2.second).Norm(), 1e-12 );
                BOOST_CHECK_LT( (box_1.second - box_3.second).Norm(), 1e-12 );

                auto delta_box = box_1.second - box_1.first;
                BOOST_CHECK_LT( std::abs(delta[0]/number_of_elements[0] - delta_box[0]) / delta_box[0], 1e-12);
                BOOST_CHECK_LT( std::abs(delta[1]/number_of_elements[1] - delta_box[1]) / delta_box[1], 1e-12);
                BOOST_CHECK_LT( std::abs(delta[2]/number_of_elements[2] - delta_box[2]) / delta_box[2], 1e-12);

                volume += delta_box[0]*delta_box[1]*delta_box[2];

                IndexType index_1 = Mapping::GetVectorIndexFromMatrixIndices(i, j, k, number_of_elements);
                IndexType index_2 = Mapping::GetVectorIndexFromMatrixIndices(indices, number_of_elements);
                BOOST_CHECK_EQUAL(index_1, index);
                BOOST_CHECK_EQUAL(index_2, index);

                PointType point_global(i+100.0, j-44.27, k+88.90);
                PointType point_param =  Mapping::PointFromGlobalToParam(point_global, lower_bound, upper_bound);
                PointType point_global_2 =  Mapping::PointFromParamToGlobal(point_param, lower_bound, upper_bound);

                BOOST_CHECK_LT( (point_global - point_global_2).Norm(), 1e-12 );
            }
        }
    }

    double volume_ref = delta[0]*delta[1]*delta[2];
    BOOST_CHECK_LT( std::abs(volume - volume_ref) / volume_ref , 1e-12);
}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra