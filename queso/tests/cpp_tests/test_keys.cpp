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
QuESo_REGISTER_KEY_SET_1( TestKeys1, SubDict, QuESo_LIST(zero, one, two, three, four) );
QuESo_REGISTER_KEY_SET_1( TestKeys2, List, QuESo_LIST(zero, one, two, three) );
QuESo_REGISTER_KEY_SET_1( TestKeys3, DataSet, QuESo_LIST(zero, one) );

BOOST_AUTO_TEST_CASE(TestRegisterKeys1) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 1" << std::endl;

    /// TestKeys1
    QuESo_CHECK_EQUAL(TestKeys1::zero, 0);
    QuESo_CHECK_EQUAL(TestKeys1::one, 1);    
    QuESo_CHECK_EQUAL(TestKeys1::two, 2);  
    QuESo_CHECK_EQUAL(TestKeys1::three, 3);
    QuESo_CHECK_EQUAL(TestKeys1::four, 4);    

    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys1::zero), "zero");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys1::one), "one");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys1::two), "two");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys1::three), "three");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys1::four), "four");
    
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys1SubDict>("zero"), TestKeys1::zero);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys1SubDict>("one"), TestKeys1::one);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys1SubDict>("two"), TestKeys1::two);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys1SubDict>("three"), TestKeys1::three);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys1SubDict>("four"), TestKeys1::four);
    BOOST_REQUIRE_THROW( key::StringToKey<key::TestKeys1SubDict>("five"), std::exception);

    typedef decltype(key::GetKeyBaseType<TestKeys1::KeyToSubDict>()) BaseType;
    static_assert( std::is_same<typename BaseType::KeyToWhat, queso::key::KeyToSubDict>::value );
    Unique<queso::key::KeyInformation> p_key_info_1 = MakeUnique<typename BaseType::KeyInfo>();
    
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyName(0), "zero");
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyName(1), "one");
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyName(2), "two");
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyName(3), "three");
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyName(4), "four");
    BOOST_REQUIRE_THROW( p_key_info_1->GetKeyName(5), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyValue("zero"), 0);
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyValue("one"), 1);
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyValue("two"), 2);
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyValue("three"), 3);
    QuESo_CHECK_EQUAL( p_key_info_1->GetKeyValue("four"), 4);
    BOOST_REQUIRE_THROW( p_key_info_1->GetKeyValue("five"), std::exception );    

    QuESo_CHECK_EQUAL( p_key_info_1->GetNumberOfKeys(), 5);

    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKeys1::zero) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKeys1::one) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKeys1::two) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKeys1::three) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, TestKeys1::four) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKeys1SubDict>("zero")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKeys1SubDict>("one")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKeys1SubDict>("two")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKeys1SubDict>("three")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_1, key::StringToKey<key::TestKeys1SubDict>("four")) );

    QuESo_CHECK( !key::IsCorrectType(p_key_info_1, TestKeys2::zero ));

    QuESo_CHECK_EQUAL( p_key_info_1->GetAllKeyNames(), "['zero', 'one', 'two', 'three', 'four']" )
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys2) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 2" << std::endl;

    /// TestKeys2
    QuESo_CHECK_EQUAL(TestKeys2::zero, 0);
    QuESo_CHECK_EQUAL(TestKeys2::one, 1);
    QuESo_CHECK_EQUAL(TestKeys2::two, 2);
    QuESo_CHECK_EQUAL(TestKeys2::three, 3);

    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys2::zero), "zero");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys2::one), "one");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys2::two), "two");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys2::three), "three");

    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys2List>("zero"), TestKeys2::zero);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys2List>("one"), TestKeys2::one);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys2List>("two"), TestKeys2::two);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys2List>("three"), TestKeys2::three);
    BOOST_REQUIRE_THROW( key::StringToKey<key::TestKeys2List>("four"), std::exception);

    typedef decltype(key::GetKeyBaseType<TestKeys2::KeyToList>()) BaseType;
    static_assert( std::is_same<typename BaseType::KeyToWhat, queso::key::KeyToList>::value );
    Unique<queso::key::KeyInformation> p_key_info_2 = MakeUnique<typename BaseType::KeyInfo>();

    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(0), "zero");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(1), "one");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(2), "two");
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyName(3), "three");
    BOOST_REQUIRE_THROW( p_key_info_2->GetKeyName(4), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("zero"), 0);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("one"), 1);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("two"), 2);
    QuESo_CHECK_EQUAL( p_key_info_2->GetKeyValue("three"), 3);
    BOOST_REQUIRE_THROW( p_key_info_2->GetKeyValue("four"), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_2->GetNumberOfKeys(), 4);

    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKeys2::zero) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKeys2::one) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKeys2::two) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, TestKeys2::three) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKeys2List>("zero")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKeys2List>("one")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKeys2List>("two")) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_2, key::StringToKey<key::TestKeys2List>("three")) );

    QuESo_CHECK( !key::IsCorrectType(p_key_info_2, TestKeys1::zero ));

    QuESo_CHECK_EQUAL( p_key_info_2->GetAllKeyNames(), "['zero', 'one', 'two', 'three']" )
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
    static_assert( std::is_same<typename BaseType::KeyToWhat, queso::key::KeyToDataSet>::value );
    Unique<queso::key::KeyInformation> p_key_info_3 = MakeUnique<typename BaseType::KeyInfo>();

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

QuESo_REGISTER_KEY_SET_2( TestKeys4, SubDict, QuESo_LIST(zero, one, two, three, four),
                                     List, QuESo_LIST(five, six, seven) );

BOOST_AUTO_TEST_CASE(TestRegisterKeys4) {
    QuESo_INFO << "Testing :: Test Keys :: Register Double Key Set" << std::endl;

    /// TestKeys4 - SubDict
    QuESo_CHECK_EQUAL(TestKeys4::zero, 0);
    QuESo_CHECK_EQUAL(TestKeys4::one, 1);
    QuESo_CHECK_EQUAL(TestKeys4::two, 2);
    QuESo_CHECK_EQUAL(TestKeys4::three, 3);
    QuESo_CHECK_EQUAL(TestKeys4::four, 4);

    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys4::zero), "zero");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys4::one), "one");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys4::two), "two");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys4::three), "three");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys4::four), "four");

    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys4SubDict>("zero"), TestKeys4::zero);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys4SubDict>("one"), TestKeys4::one);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys4SubDict>("two"), TestKeys4::two);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys4SubDict>("three"), TestKeys4::three);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys4SubDict>("four"), TestKeys4::four);
    BOOST_REQUIRE_THROW( key::StringToKey<key::TestKeys4SubDict>("five"), std::exception);

    /// TestKeys4 - List
    QuESo_CHECK_EQUAL(TestKeys4::five, 0);
    QuESo_CHECK_EQUAL(TestKeys4::six, 1);
    QuESo_CHECK_EQUAL(TestKeys4::seven, 2);

    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys4::five), "five");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys4::six), "six");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys4::seven), "seven");

    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys4List>("five"), TestKeys4::five);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys4List>("six"), TestKeys4::six);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys4List>("seven"), TestKeys4::seven);
    BOOST_REQUIRE_THROW( key::StringToKey<key::TestKeys4List>("eight"), std::exception);

    typedef decltype(key::GetKeyBaseType<TestKeys4::KeyToSubDict>()) BaseTypeSubDict;
    static_assert( std::is_same<typename BaseTypeSubDict::KeyToWhat, queso::key::KeyToSubDict>::value );
    Unique<queso::key::KeyInformation> p_key_info_4_subdict = MakeUnique<typename BaseTypeSubDict::KeyInfo>();

    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetKeyName(0), "zero");
    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetKeyName(1), "one");
    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetKeyName(2), "two");
    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetKeyName(3), "three");
    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetKeyName(4), "four");
    BOOST_REQUIRE_THROW( p_key_info_4_subdict->GetKeyName(5), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetKeyValue("zero"), 0);
    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetKeyValue("one"), 1);
    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetKeyValue("two"), 2);
    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetKeyValue("three"), 3);
    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetKeyValue("four"), 4);
    BOOST_REQUIRE_THROW( p_key_info_4_subdict->GetKeyValue("five"), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetNumberOfKeys(), 5);

    typedef decltype(key::GetKeyBaseType<TestKeys4::KeyToList>()) BaseTypeList;
    static_assert( std::is_same<typename BaseTypeList::KeyToWhat, queso::key::KeyToList>::value );
    Unique<queso::key::KeyInformation> p_key_info_4_list = MakeUnique<typename BaseTypeList::KeyInfo>();

    QuESo_CHECK_EQUAL( p_key_info_4_list->GetKeyName(0), "five");
    QuESo_CHECK_EQUAL( p_key_info_4_list->GetKeyName(1), "six");
    QuESo_CHECK_EQUAL( p_key_info_4_list->GetKeyName(2), "seven");
    BOOST_REQUIRE_THROW( p_key_info_4_list->GetKeyName(3), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_4_list->GetKeyValue("five"), 0);
    QuESo_CHECK_EQUAL( p_key_info_4_list->GetKeyValue("six"), 1);
    QuESo_CHECK_EQUAL( p_key_info_4_list->GetKeyValue("seven"), 2);
    BOOST_REQUIRE_THROW( p_key_info_4_list->GetKeyValue("eight"), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_4_list->GetNumberOfKeys(), 3);

    QuESo_CHECK( key::IsCorrectType(p_key_info_4_subdict, TestKeys4::zero) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_4_subdict, TestKeys4::one) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_4_subdict, TestKeys4::two) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_4_subdict, TestKeys4::three) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_4_subdict, TestKeys4::four) );

    QuESo_CHECK( key::IsCorrectType(p_key_info_4_list, TestKeys4::five) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_4_list, TestKeys4::six) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_4_list, TestKeys4::seven) );

    QuESo_CHECK( !key::IsCorrectType(p_key_info_4_subdict, TestKeys2::zero) );
    QuESo_CHECK( !key::IsCorrectType(p_key_info_4_list, TestKeys1::zero) );

    QuESo_CHECK_EQUAL( p_key_info_4_subdict->GetAllKeyNames(), "['zero', 'one', 'two', 'three', 'four']" );
    QuESo_CHECK_EQUAL( p_key_info_4_list->GetAllKeyNames(), "['five', 'six', 'seven']" );
}

QuESo_REGISTER_KEY_SET_3( TestKeys5, SubDict, QuESo_LIST(one, two, three, four),
                                     List, QuESo_LIST(five, six, seven),
                                     DataSet, QuESo_LIST(eight, nine, ten) );

BOOST_AUTO_TEST_CASE(TestRegisterKeys5) {
    QuESo_INFO << "Testing :: Test Keys :: Register Triple Key Set" << std::endl;

    /// TestKeys5 - SubDict
    QuESo_CHECK_EQUAL(TestKeys5::one, 0);
    QuESo_CHECK_EQUAL(TestKeys5::two, 1);
    QuESo_CHECK_EQUAL(TestKeys5::three, 2);
    QuESo_CHECK_EQUAL(TestKeys5::four, 3);

    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys5::one), "one");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys5::two), "two");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys5::three), "three");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys5::four), "four");

    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys5SubDict>("one"), TestKeys5::one);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys5SubDict>("two"), TestKeys5::two);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys5SubDict>("three"), TestKeys5::three);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys5SubDict>("four"), TestKeys5::four);
    BOOST_REQUIRE_THROW( key::StringToKey<key::TestKeys5SubDict>("five"), std::exception);
    /// TestKeys5 - List
    QuESo_CHECK_EQUAL(TestKeys5::five, 0);
    QuESo_CHECK_EQUAL(TestKeys5::six, 1);
    QuESo_CHECK_EQUAL(TestKeys5::seven, 2);

    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys5::five), "five");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys5::six), "six");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys5::seven), "seven");

    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys5List>("five"), TestKeys5::five);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys5List>("six"), TestKeys5::six);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys5List>("seven"), TestKeys5::seven);
    BOOST_REQUIRE_THROW( key::StringToKey<key::TestKeys5List>("eight"), std::exception);

    /// TestKeys5 - DataSet
    QuESo_CHECK_EQUAL(TestKeys5::eight, 0);
    QuESo_CHECK_EQUAL(TestKeys5::nine, 1);
    QuESo_CHECK_EQUAL(TestKeys5::ten, 2);

    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys5::eight), "eight");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys5::nine), "nine");
    QuESo_CHECK_EQUAL( key::KeyToString(TestKeys5::ten), "ten");

    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys5DataSet>("eight"), TestKeys5::eight);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys5DataSet>("nine"), TestKeys5::nine);
    QuESo_CHECK_EQUAL( key::StringToKey<key::TestKeys5DataSet>("ten"), TestKeys5::ten);
    BOOST_REQUIRE_THROW( key::StringToKey<key::TestKeys5DataSet>("eleven"), std::exception);

    typedef decltype(key::GetKeyBaseType<TestKeys5::KeyToSubDict>()) BaseTypeSubDict;
    static_assert( std::is_same<typename BaseTypeSubDict::KeyToWhat, queso::key::KeyToSubDict>::value );
    Unique<queso::key::KeyInformation> p_key_info_5_subdict = MakeUnique<typename BaseTypeSubDict::KeyInfo>();

    QuESo_CHECK_EQUAL( p_key_info_5_subdict->GetKeyName(0), "one");
    QuESo_CHECK_EQUAL( p_key_info_5_subdict->GetKeyName(1), "two");
    QuESo_CHECK_EQUAL( p_key_info_5_subdict->GetKeyName(2), "three");
    QuESo_CHECK_EQUAL( p_key_info_5_subdict->GetKeyName(3), "four");
    BOOST_REQUIRE_THROW( p_key_info_5_subdict->GetKeyName(4), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_5_subdict->GetKeyValue("one"), 0);
    QuESo_CHECK_EQUAL( p_key_info_5_subdict->GetKeyValue("two"), 1);
    QuESo_CHECK_EQUAL( p_key_info_5_subdict->GetKeyValue("three"), 2);
    QuESo_CHECK_EQUAL( p_key_info_5_subdict->GetKeyValue("four"), 3);
    BOOST_REQUIRE_THROW( p_key_info_5_subdict->GetKeyValue("five"), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_5_subdict->GetNumberOfKeys(), 4);

    typedef decltype(key::GetKeyBaseType<TestKeys5::KeyToList>()) BaseTypeList;
    static_assert( std::is_same<typename BaseTypeList::KeyToWhat, queso::key::KeyToList>::value );
    Unique<queso::key::KeyInformation> p_key_info_5_list = MakeUnique<typename BaseTypeList::KeyInfo>();

    QuESo_CHECK_EQUAL( p_key_info_5_list->GetKeyName(0), "five");
    QuESo_CHECK_EQUAL( p_key_info_5_list->GetKeyName(1), "six");
    QuESo_CHECK_EQUAL( p_key_info_5_list->GetKeyName(2), "seven");
    BOOST_REQUIRE_THROW( p_key_info_5_list->GetKeyName(3), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_5_list->GetKeyValue("five"), 0);
    QuESo_CHECK_EQUAL( p_key_info_5_list->GetKeyValue("six"), 1);
    QuESo_CHECK_EQUAL( p_key_info_5_list->GetKeyValue("seven"), 2);
    BOOST_REQUIRE_THROW( p_key_info_5_list->GetKeyValue("eight"), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_5_list->GetNumberOfKeys(), 3);

    typedef decltype(key::GetKeyBaseType<TestKeys5::KeyToDataSet>()) BaseTypeDataSet;
    static_assert( std::is_same<typename BaseTypeDataSet::KeyToWhat, queso::key::KeyToDataSet>::value );
    Unique<queso::key::KeyInformation> p_key_info_5_dataset = MakeUnique<typename BaseTypeDataSet::KeyInfo>();

    QuESo_CHECK_EQUAL( p_key_info_5_dataset->GetKeyName(0), "eight");
    QuESo_CHECK_EQUAL( p_key_info_5_dataset->GetKeyName(1), "nine");
    QuESo_CHECK_EQUAL( p_key_info_5_dataset->GetKeyName(2), "ten");
    BOOST_REQUIRE_THROW( p_key_info_5_dataset->GetKeyName(3), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_5_dataset->GetKeyValue("eight"), 0);
    QuESo_CHECK_EQUAL( p_key_info_5_dataset->GetKeyValue("nine"), 1);
    QuESo_CHECK_EQUAL( p_key_info_5_dataset->GetKeyValue("ten"), 2);
    BOOST_REQUIRE_THROW( p_key_info_5_dataset->GetKeyValue("eleven"), std::exception );

    QuESo_CHECK_EQUAL( p_key_info_5_dataset->GetNumberOfKeys(), 3);

    QuESo_CHECK( key::IsCorrectType(p_key_info_5_subdict, TestKeys5::one) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_5_subdict, TestKeys5::two) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_5_subdict, TestKeys5::three) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_5_subdict, TestKeys5::four) );

    QuESo_CHECK( key::IsCorrectType(p_key_info_5_list, TestKeys5::five) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_5_list, TestKeys5::six) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_5_list, TestKeys5::seven) );

    QuESo_CHECK( key::IsCorrectType(p_key_info_5_dataset, TestKeys5::eight) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_5_dataset, TestKeys5::nine) );
    QuESo_CHECK( key::IsCorrectType(p_key_info_5_dataset, TestKeys5::ten) );

    QuESo_CHECK( !key::IsCorrectType(p_key_info_5_subdict, TestKeys2::zero) );
    QuESo_CHECK( !key::IsCorrectType(p_key_info_5_list, TestKeys1::zero) );
    QuESo_CHECK( !key::IsCorrectType(p_key_info_5_dataset, TestKeys3::zero) );

    QuESo_CHECK_EQUAL( p_key_info_5_subdict->GetAllKeyNames(), "['one', 'two', 'three', 'four']" );
    QuESo_CHECK_EQUAL( p_key_info_5_list->GetAllKeyNames(), "['five', 'six', 'seven']" );
    QuESo_CHECK_EQUAL( p_key_info_5_dataset->GetAllKeyNames(), "['eight', 'nine', 'ten']" );
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso