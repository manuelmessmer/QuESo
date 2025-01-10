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

// Due to the Testing namespace, we have to establish the following typedefs.
namespace key {
typedef queso::key::SubDict SubDict;
typedef queso::key::List List;
typedef queso::key::DataSet DataSet;
typedef queso::key::KeyInformation KeyInformation;
}

QuESo_REGISTER_KEY_SET_1( TestKeys1, SubDict, QuESo_LIST(zero, one, three=3, four, six=6) );
QuESo_REGISTER_KEY_SET_1( TestKeys2, List, QuESo_LIST(zero, four=4, five, six, seven) );
QuESo_REGISTER_KEY_SET_1( TestKeys3, DataSet, QuESo_LIST(zero, one) );

BOOST_AUTO_TEST_CASE(TestRegisterKeys1) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 1" << std::endl;

    /// TestKeys1
    QuESo_CHECK_EQUAL(TestKeys1::zero, 0);
    QuESo_CHECK_EQUAL(TestKeys1::one, 1);    
    QuESo_CHECK_EQUAL(TestKeys1::three, 3);  
    QuESo_CHECK_EQUAL(TestKeys1::four, 4);
    QuESo_CHECK_EQUAL(TestKeys1::six, 6);    

    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys1::zero), "zero");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys1::one), "one");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys1::three), "three");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys1::four), "four");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys1::six), "six");
    
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys1SubDict>("zero"), TestKeys1::zero);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys1SubDict>("one"), TestKeys1::one);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys1SubDict>("three"), TestKeys1::three);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys1SubDict>("four"), TestKeys1::four);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys1SubDict>("six"), TestKeys1::six);
    BOOST_REQUIRE_THROW( key::StringToKey<key::TestKeys1SubDict>("two"), std::exception);

    typedef decltype(key::GetKeyBaseType<TestKeys1::KeyToSubDict>()) BaseType;
    static_assert( std::is_same<typename BaseType::KeyToWhat, key::SubDict>::value );
    Unique<key::KeyInformation> p_key_info_1 = MakeUnique<typename BaseType::KeyInfo>();
    
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

    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKeys1::zero) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKeys1::one) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKeys1::three) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKeys1::four) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKeys1::six) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKeys1SubDict>("zero")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKeys1SubDict>("one")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKeys1SubDict>("three")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKeys1SubDict>("four")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKeys1SubDict>("six")) );

    QuESo_CHECK( !key::IsCorrectType(p_key_info_1, TestKeys2::zero ));

    QuESo_CHECK_EQUAL( p_key_info_1->GetAllKeyNames(), "['zero', 'one', 'three', 'four', 'six']" )
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys2) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 2" << std::endl;

    /// TestKeys2
    QuESo_CHECK_EQUAL(TestKeys2::zero, 0);
    QuESo_CHECK_EQUAL(TestKeys2::four, 4);    
    QuESo_CHECK_EQUAL(TestKeys2::five, 5);  
    QuESo_CHECK_EQUAL(TestKeys2::six, 6);
    QuESo_CHECK_EQUAL(TestKeys2::seven, 7);    

    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys2::zero), "zero");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys2::four), "four");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys2::five), "five");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys2::six), "six");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys2::seven), "seven");
    
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys2List>("zero"), TestKeys2::zero);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys2List>("four"), TestKeys2::four);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys2List>("five"), TestKeys2::five);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys2List>("six"), TestKeys2::six);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys2List>("seven"), TestKeys2::seven);
    BOOST_REQUIRE_THROW( key::StringToKey<key::TestKeys2List>("three"), std::exception);

    typedef decltype(key::GetKeyBaseType<TestKeys2::KeyToList>()) BaseType;
    static_assert( std::is_same<typename BaseType::KeyToWhat, key::List>::value );
    Unique<key::KeyInformation> p_key_info_2 = MakeUnique<typename BaseType::KeyInfo>();
    
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(0), "zero");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(4), "four");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(5), "five");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(6), "six");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(7), "seven");
    BOOST_REQUIRE_THROW( p_key_info_2->GetKeyName(3), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("zero"), 0);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("four"), 4);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("five"), 5);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("six"), 6);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("seven"), 7);
    BOOST_REQUIRE_THROW( p_key_info_2->GetKeyValue("three"), std::exception );    

    QuESo_CHECK_EQUAL( p_key_info_2->GetNumberOfKeys(), 5);

    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKeys2::zero) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKeys2::four) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKeys2::five) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKeys2::six) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKeys2::seven) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKeys2List>("zero")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKeys2List>("four")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKeys2List>("five")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKeys2List>("six")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKeys2List>("seven")) );

    QuESo_CHECK( !key::IsCorrectType(p_key_info_2, TestKeys1::zero ));

    QuESo_CHECK_EQUAL( p_key_info_2->GetAllKeyNames(), "['zero', 'four', 'five', 'six', 'seven']" )
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys3) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 3" << std::endl;

    /// TestKeys3
    QuESo_CHECK_EQUAL(TestKeys3::zero, 0);
    QuESo_CHECK_EQUAL(TestKeys3::one, 1);

    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys3::zero), "zero");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys3::one), "one");

    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys3DataSet>("zero"), TestKeys3::zero);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys3DataSet>("one"), TestKeys3::one);
    BOOST_REQUIRE_THROW( key::StringToKey<key::TestKeys3DataSet>("two"), std::exception);

    typedef decltype(key::GetKeyBaseType<TestKeys3::KeyToDataSet>()) BaseType;
    static_assert( std::is_same<typename BaseType::KeyToWhat, key::DataSet>::value );
    Unique<key::KeyInformation> p_key_info_3 = MakeUnique<typename BaseType::KeyInfo>();

    QuESo_CHECK_EQUAL( p_key_info_3->GetKeyName(0), "zero");
    QuESo_CHECK_EQUAL( p_key_info_3->GetKeyName(1), "one");
    BOOST_REQUIRE_THROW( p_key_info_3->GetKeyName(2), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_3->GetKeyValue("zero"), 0);
    QuESo_CHECK_EQUAL( p_key_info_3->GetKeyValue("one"), 1);
    BOOST_REQUIRE_THROW( p_key_info_3->GetKeyValue("two"), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_3->GetNumberOfKeys(), 2);

    QuESo_CHECK( key::IsCorrectType(p_key_info_3, TestKeys3::zero) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_3, TestKeys3::one) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_3, key::StringToKey<key::TestKeys3DataSet>("zero")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_3, key::StringToKey<key::TestKeys3DataSet>("one")) );

    QuESo_CHECK( !key::IsCorrectType(p_key_info_3, TestKeys1::zero) );
    QuESo_CHECK( !key::IsCorrectType(p_key_info_3, TestKeys2::zero) );

    QuESo_CHECK_EQUAL( p_key_info_3->GetAllKeyNames(), "['zero', 'one']" )
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso