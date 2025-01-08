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
namespace Testing {

BOOST_AUTO_TEST_SUITE( KeysTestSuite )

// Due to th Testing namespace, we have to establish the following typedefs.
namespace key {
typedef queso::key::SubDict SubDict;
typedef queso::key::List List;
typedef queso::key::DataSet DataSet;
typedef queso::key::KeyInformation KeyInformation;
}

QuESo_REGISTER_KEYS_1( TestKey1, SubDict, QuESo_LIST(zero, one, three=3, four, six=6) );
QuESo_REGISTER_KEYS_1( TestKey2, List, QuESo_LIST(zero, four=4, five, six, seven) );
QuESo_REGISTER_KEYS_1( TestKey3, DataSet, QuESo_LIST(zero, one) );

BOOST_AUTO_TEST_CASE(TestRegisterKeys1) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 1" << std::endl;

    /// TestKey1
    QuESo_CHECK_EQUAL(TestKey1::zero, 0);
    QuESo_CHECK_EQUAL(TestKey1::one, 1);    
    QuESo_CHECK_EQUAL(TestKey1::three, 3);  
    QuESo_CHECK_EQUAL(TestKey1::four, 4);
    QuESo_CHECK_EQUAL(TestKey1::six, 6);    

    QuESo_CHECK_EQUAL( key::KeyToString(TestKey1::zero), "zero");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKey1::one), "one");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKey1::three), "three");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKey1::four), "four");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKey1::six), "six");
    
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKey1SubDict>("zero"), TestKey1::zero);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKey1SubDict>("one"), TestKey1::one);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKey1SubDict>("three"), TestKey1::three);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKey1SubDict>("four"), TestKey1::four);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKey1SubDict>("six"), TestKey1::six);

    typedef decltype(key::GetKeyBaseType<TestKey1::KeyToSubDict>()) BaseType1;
    static_assert( std::is_same<typename BaseType1::KeyToWhat, key::SubDict>::value );
    Unique<key::KeyInformation> p_key_info_1 = MakeUnique<typename BaseType1::KeyInfo>();
    
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyName(0), "zero");
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyName(1), "one");
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyName(3), "three");
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyName(4), "four");
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyName(6), "six");
    BOOST_REQUIRE_THROW( p_key_info_1->GetKeyName(2), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyValue("zero"), 0);
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyValue("one"), 1);
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyValue("three"), 3);
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyValue("four"), 4);
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyValue("six"), 6);
    BOOST_REQUIRE_THROW( p_key_info_1->GetKeyValue("two"), std::exception );    

    QuESo_CHECK_EQUAL( p_key_info_1->GetNumberOfKeys(), 5);

    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKey1::zero) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKey1::one) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKey1::three) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKey1::four) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKey1::six) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKey1SubDict>("zero")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKey1SubDict>("one")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKey1SubDict>("three")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKey1SubDict>("four")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKey1SubDict>("six")) );

    QuESo_CHECK( !key::IsCorrectType(p_key_info_1, TestKey2::zero ));

    QuESo_CHECK_EQUAL( p_key_info_1->GetAllKeyNames(), "[zero, one, three, four, six]" )
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys2) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 2" << std::endl;

    /// TestKey2
    QuESo_CHECK_EQUAL(TestKey2::zero, 0);
    QuESo_CHECK_EQUAL(TestKey2::four, 4);
    QuESo_CHECK_EQUAL(TestKey2::five, 5);
    QuESo_CHECK_EQUAL(TestKey2::six, 6);
    QuESo_CHECK_EQUAL(TestKey2::seven, 7);

    QuESo_CHECK_EQUAL( key::KeyToString(TestKey2::zero), "zero");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKey2::four), "four");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKey2::five), "five");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKey2::six), "six");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKey2::seven), "seven");

    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKey2List>("zero"), TestKey2::zero);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKey2List>("four"), TestKey2::four);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKey2List>("five"), TestKey2::five);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKey2List>("six"), TestKey2::six);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKey2List>("seven"), TestKey2::seven);

    typedef decltype(key::GetKeyBaseType<TestKey2::KeyToList>()) BaseType2;
    static_assert( std::is_same<typename BaseType2::KeyToWhat, key::List>::value );
    Unique<key::KeyInformation> p_key_info_2 = MakeUnique<typename BaseType2::KeyInfo>();

    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(0), "zero");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(4), "four");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(5), "five");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(6), "six");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(7), "seven");
    BOOST_REQUIRE_THROW( p_key_info_2->GetKeyName(1), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("zero"), 0);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("four"), 4);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("five"), 5);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("six"), 6);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("seven"), 7);
    BOOST_REQUIRE_THROW( p_key_info_2->GetKeyValue("one"), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_2->GetNumberOfKeys(), 5);

    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKey2::zero) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKey2::four) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKey2::five) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKey2::six) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKey2::seven) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKey2List>("zero")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKey2List>("four")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKey2List>("five")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKey2List>("six")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKey2List>("seven")) );

    QuESo_CHECK( !key::IsCorrectType(p_key_info_2, TestKey1::zero ));

    QuESo_CHECK_EQUAL( p_key_info_2->GetAllKeyNames(), "[zero, four, five, six, seven]" )
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso