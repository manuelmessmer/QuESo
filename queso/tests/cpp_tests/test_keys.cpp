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


namespace queso {

QuESo_REGISTER_KEYS_1( TestKey1, SubDict, QuESo_LIST(zero, one, three=3, four, six=6) );

namespace Testing {
    

BOOST_AUTO_TEST_SUITE( KeysTestSuite )

BOOST_AUTO_TEST_CASE(KeysMacros) {
    QuESo_INFO << "Testing :: Test keys :: Macros" << std::endl;

    auto peter_enum_one = TestKey1::one;
    static_assert( std::is_same<decltype(GetType(peter_enum_one)), key::SubDict>::value );

    QuESo_CHECK_NOT_EQUAL(TestKey1::one, TestKey1::zero);

    QuESo_CHECK_EQUAL(key::KeyToString(TestKey1::zero), "zero");
    QuESo_CHECK_EQUAL(key::KeyToString(TestKey1::one), "one");
    
    // auto peter_enum_two = key::StringToKey<Peter::BaseEnum>("two");

    TestKey1 peter_enum_info{};
    key::KeyInformation* p_key_info = &peter_enum_info;
    QuESo_CHECK_EQUAL( p_key_info->EnumName(0), "zero");
    QuESo_CHECK_EQUAL( p_key_info->EnumName(1), "one");
    QuESo_CHECK_EQUAL( p_key_info->EnumName(3), "three");
    QuESo_CHECK_EQUAL( p_key_info->EnumName(4), "four");
    QuESo_CHECK_EQUAL( p_key_info->EnumName(6), "six");
    
}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso