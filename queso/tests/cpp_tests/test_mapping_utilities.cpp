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
//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/containers/grid_indexer.hpp"
#include "queso/utilities/mapping_utilities.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( MappingTestSuite )

BOOST_AUTO_TEST_CASE(MappingTest) {
    QuESo_INFO << "Testing :: Test Mapping Utilities :: Mapping" << std::endl;

    const BoundingBoxType bounds_xyz = MakeBox( {-25, -111.44, 7.89}, {78.67, -35.68, 18.99} );
    const BoundingBoxType bounds_uvw = MakeBox( { -10.0,  -2.2, 2.0}, {2.0, 10.0, 17.0} );

    const Vector3i number_of_elements{15, 11, 17};

    for( IndexType i = 0; i < number_of_elements[0]; ++i){
        for( IndexType j = 0; j < number_of_elements[1]; ++j){
            for( IndexType k = 0; k < number_of_elements[2]; ++k){

                PointType point_global{i+100.0, j-44.27, k+88.90};
                PointType point_param =  Mapping::PointFromGlobalToParam(point_global, bounds_xyz, bounds_uvw);
                PointType point_global_2 =  Mapping::PointFromParamToGlobal(point_param, bounds_xyz, bounds_uvw);

                QuESo_CHECK_LT( Math::Norm( Math::Subtract(point_global, point_global_2) ), 1e-12 );
            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso
