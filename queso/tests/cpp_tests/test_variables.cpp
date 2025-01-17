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
#include "queso/includes/keys.hpp"
#include "queso/includes/variables.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( VariablesTestSuite )

// Due to the Testing namespace, we have to establish the following typedefs.
namespace key {
typedef queso::key::SubDict SubDict;
typedef queso::key::List List;
typedef queso::key::DataSet DataSet;
typedef queso::key::KeyInformation KeyInformation;
}
namespace variable {
namespace detail {
using queso::variable::detail::CheckIfAllKeysAreSet;
}
}

QuESo_REGISTER_KEY_SET_1( TestKeys1, DataSet, 
    QuESo_LIST(zero, one, two, three, four, five, six, seven) );

QuESo_REGISTER_DATASET_VARIABLES(TestKeys1, 
    QuESo_VARIABLE(TestKeys1::zero, PointType),
    QuESo_VARIABLE(TestKeys1::one, Vector3i),
    QuESo_VARIABLE(TestKeys1::two, bool),
    QuESo_VARIABLE(TestKeys1::three, double),
    QuESo_VARIABLE(TestKeys1::four, IndexType),
    QuESo_VARIABLE(TestKeys1::five, std::string),
    QuESo_VARIABLE(TestKeys1::six, IntegrationMethodType),
    QuESo_VARIABLE(TestKeys1::seven, GridTypeType),
);

BOOST_AUTO_TEST_CASE(TestVariablesValueTypes1) {
    QuESo_INFO << "Testing :: Test Variables :: Chech Value Types Single Set" << std::endl;
  
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys1::zero)) == std::type_index(typeid(PointType)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys1::one)) == std::type_index(typeid(Vector3i)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys1::two)) == std::type_index(typeid(bool)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys1::three)) == std::type_index(typeid(double)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys1::four)) == std::type_index(typeid(IndexType)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys1::five)) == std::type_index(typeid(std::string)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys1::six)) == std::type_index(typeid(IntegrationMethodType)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys1::seven)) == std::type_index(typeid(GridTypeType)) );
}

QuESo_REGISTER_KEY_SET_2( TestKeys2, SubDict, QuESo_LIST(zero, one),
                                     DataSet, QuESo_LIST(two, three, four, five, six, seven, eight, nine) );

QuESo_REGISTER_DATASET_VARIABLES(TestKeys2, 
    QuESo_VARIABLE(TestKeys2::two, PointType),
    QuESo_VARIABLE(TestKeys2::three, Vector3i),
    QuESo_VARIABLE(TestKeys2::four, bool),
    QuESo_VARIABLE(TestKeys2::five, double),
    QuESo_VARIABLE(TestKeys2::six, IndexType),
    QuESo_VARIABLE(TestKeys2::seven, std::string),
    QuESo_VARIABLE(TestKeys2::eight, IntegrationMethodType),
    QuESo_VARIABLE(TestKeys2::nine, GridTypeType),
);         

BOOST_AUTO_TEST_CASE(TestVariablesValueTypes2) {
    QuESo_INFO << "Testing :: Test Variables :: Check Value Types Double Sets" << std::endl;
    
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys2::two)) == std::type_index(typeid(PointType)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys2::three)) == std::type_index(typeid(Vector3i)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys2::four)) == std::type_index(typeid(bool)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys2::five)) == std::type_index(typeid(double)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys2::six)) == std::type_index(typeid(IndexType)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys2::seven)) == std::type_index(typeid(std::string)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys2::eight)) == std::type_index(typeid(IntegrationMethodType)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys2::nine)) == std::type_index(typeid(GridTypeType)) );
}

QuESo_REGISTER_KEY_SET_3( TestKeys3, SubDict, QuESo_LIST(zero, one),
                                     DataSet, QuESo_LIST(two, three, four, five, six, seven, eight, nine),
                                     List, QuESo_LIST(ten, eleven) );

QuESo_REGISTER_DATASET_VARIABLES(TestKeys3, 
    QuESo_VARIABLE(TestKeys3::two, PointType),
    QuESo_VARIABLE(TestKeys3::three, Vector3i),
    QuESo_VARIABLE(TestKeys3::four, bool),
    QuESo_VARIABLE(TestKeys3::five, double),
    QuESo_VARIABLE(TestKeys3::six, IndexType),
    QuESo_VARIABLE(TestKeys3::seven, std::string),
    QuESo_VARIABLE(TestKeys3::eight, IntegrationMethodType),
    QuESo_VARIABLE(TestKeys3::nine, GridTypeType),
);

BOOST_AUTO_TEST_CASE(TestVariablesValueTypes3) {
    QuESo_INFO << "Testing :: Test Variables :: Check Value Types Triple Set" << std::endl;
    
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys3::two)) == std::type_index(typeid(PointType)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys3::three)) == std::type_index(typeid(Vector3i)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys3::four)) == std::type_index(typeid(bool)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys3::five)) == std::type_index(typeid(double)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys3::six)) == std::type_index(typeid(IndexType)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys3::seven)) == std::type_index(typeid(std::string)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys3::eight)) == std::type_index(typeid(IntegrationMethodType)) );
    QuESo_CHECK( (variable::GetValueTypeIndex(TestKeys3::nine)) == std::type_index(typeid(GridTypeType)) );
}

QuESo_REGISTER_KEY_SET_1( TestKeys4, DataSet, 
    QuESo_LIST(zero, one, two, three, four, five, six, seven) );

QuESo_REGISTER_DATASET_VARIABLES(TestKeys4, 
    QuESo_VARIABLE(TestKeys4::zero, PointType),
    QuESo_VARIABLE(TestKeys4::one, Vector3i),
    QuESo_VARIABLE(TestKeys4::two, bool),
    QuESo_VARIABLE(TestKeys4::three, double),
    QuESo_VARIABLE(TestKeys4::four, IndexType),
    QuESo_VARIABLE(TestKeys4::five, std::string)
);

BOOST_AUTO_TEST_CASE(TestVariablesException1) {
    QuESo_INFO << "Testing :: Test Variables :: Check Exceptions 1" << std::endl;

    if( !NOTDEBUG ) {
        BOOST_REQUIRE_THROW( variable::GetValueTypeIndex(TestKeys4::zero), std::exception );
    }
}

QuESo_REGISTER_KEY_SET_1( TestKeys5, DataSet, 
    QuESo_LIST(zero, one, two, three, four, five, six, seven) );

QuESo_REGISTER_DATASET_VARIABLES(TestKeys5, 
    QuESo_VARIABLE(TestKeys5::zero, PointType),
    QuESo_VARIABLE(TestKeys5::one, Vector3i),
    QuESo_VARIABLE(TestKeys5::two, bool),
    QuESo_VARIABLE(TestKeys5::three, double),
    QuESo_VARIABLE(TestKeys5::four, IndexType),
    QuESo_VARIABLE(TestKeys5::five, std::string),
    QuESo_VARIABLE(TestKeys5::six, IntegrationMethodType),
    QuESo_VARIABLE(TestKeys5::six, GridTypeType),
);

BOOST_AUTO_TEST_CASE(TestVariablesException5) {
    QuESo_INFO << "Testing :: Test Variables :: Check Exceptions 2" << std::endl;

    if( !NOTDEBUG ) {
        BOOST_REQUIRE_THROW( variable::GetValueTypeIndex(TestKeys5::six), std::exception );
    }
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso