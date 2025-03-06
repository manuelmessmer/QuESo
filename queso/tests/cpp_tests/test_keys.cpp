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
#include "queso/includes/register_keys.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( KeysTestSuite )

template<typename TKeySetToWhat, typename TKeyToWhat = TKeySetToWhat, typename TKeyType>
void CheckKey(const TKeyType& rKey, IndexType Index, const std::string& rName) {
    QuESo_CHECK_EQUAL(rKey.Index(), Index);
    QuESo_CHECK_EQUAL(rKey.Name(), rName);
    QuESo_CHECK(rKey.VariableTypeIndex() == std::type_index(typeid(TKeyToWhat)));

    static_assert( std::is_same<typename TKeyType::KeyToWhat, TKeyToWhat>::value );
    static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeySetToWhat, TKeySetToWhat>::value );
    typename TKeyType::KeySetInfoType key_set_info;
    key_set_info.IsCorrectKeyType(rKey);

    const auto p_key = key_set_info.pGetKey(rKey.Name());
    QuESo_CHECK_EQUAL(p_key->Name(), rKey.Name());
    QuESo_CHECK_EQUAL(p_key->Index(), rKey.Index());
    QuESo_CHECK(p_key->VariableTypeIndex() == rKey.VariableTypeIndex());
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys1) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 1" << std::endl;

    /// TestKeys1
    CheckKey<queso::key::SubDictTypeTag>(TestKeys1::zero, 0, "zero");
    CheckKey<queso::key::SubDictTypeTag>(TestKeys1::one, 1, "one");
    CheckKey<queso::key::SubDictTypeTag>(TestKeys1::two, 2, "two");
    CheckKey<queso::key::SubDictTypeTag>(TestKeys1::three, 3, "three");
    CheckKey<queso::key::SubDictTypeTag>(TestKeys1::four, 4, "four");

    typedef decltype(TestKeys1::zero)::KeySetInfoType KeySetInfoType;
    KeySetInfoType key_set_info{};
    QuESo_CHECK_EQUAL( key_set_info.GetNumberOfKeys(), 5);
    QuESo_CHECK_EQUAL( KeySetInfoType::StaticGetAllKeyNames(), "['zero', 'one', 'two', 'three', 'four']" );
    BOOST_REQUIRE_THROW( key_set_info.pGetKey("five"), std::exception );
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys2) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 2" << std::endl;

    /// TestKeys2
    CheckKey<queso::key::ListTypeTag>(TestKeys2::zero, 0, "zero");
    CheckKey<queso::key::ListTypeTag>(TestKeys2::one, 1, "one");
    CheckKey<queso::key::ListTypeTag>(TestKeys2::two, 2, "two");
    CheckKey<queso::key::ListTypeTag>(TestKeys2::three, 3, "three");

    typedef decltype(TestKeys2::zero)::KeySetInfoType KeySetInfoType;
    KeySetInfoType key_set_info{};
    QuESo_CHECK_EQUAL( key_set_info.GetNumberOfKeys(), 4);
    QuESo_CHECK_EQUAL( KeySetInfoType::StaticGetAllKeyNames(), "['zero', 'one', 'two', 'three']" );
    BOOST_REQUIRE_THROW( key_set_info.pGetKey("four"), std::exception );
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys3) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 3" << std::endl;

    /// TestKeys3
    CheckKey<queso::key::MainValuesTypeTag, PointType>(TestKeys3::zero, 0, "zero");
    CheckKey<queso::key::MainValuesTypeTag, Vector3i>(TestKeys3::one, 1, "one");
    CheckKey<queso::key::MainValuesTypeTag, bool>(TestKeys3::two, 2, "two");
    CheckKey<queso::key::MainValuesTypeTag, double>(TestKeys3::three, 3, "three");
    CheckKey<queso::key::MainValuesTypeTag, IndexType>(TestKeys3::four, 4, "four");
    CheckKey<queso::key::MainValuesTypeTag, std::string>(TestKeys3::five, 5, "five");
    CheckKey<queso::key::MainValuesTypeTag, IntegrationMethodType>(TestKeys3::six, 6, "six");
    CheckKey<queso::key::MainValuesTypeTag, GridTypeType>(TestKeys3::seven, 7, "seven");

    typedef decltype(TestKeys3::zero)::KeySetInfoType KeySetInfoType;
    KeySetInfoType key_set_info{};
    QuESo_CHECK_EQUAL( key_set_info.GetNumberOfKeys(), 8);
    QuESo_CHECK_EQUAL( KeySetInfoType::StaticGetAllKeyNames(), "['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven']" );
    BOOST_REQUIRE_THROW( key_set_info.pGetKey("twelve"), std::exception );
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys4) {
    QuESo_INFO << "Testing :: Test Keys :: Register Double Key Set" << std::endl;

    /// TestKeys4 :: SubDictTypeTag
    CheckKey<queso::key::SubDictTypeTag>(TestKeys4::zero, 0, "zero");
    CheckKey<queso::key::SubDictTypeTag>(TestKeys4::one, 1, "one");
    CheckKey<queso::key::SubDictTypeTag>(TestKeys4::two, 2, "two");
    CheckKey<queso::key::SubDictTypeTag>(TestKeys4::three, 3, "three");

    typedef decltype(TestKeys4::zero)::KeySetInfoType KeySetInfoTypeSubDict;
    KeySetInfoTypeSubDict key_set_info_subdict{};
    QuESo_CHECK_EQUAL( key_set_info_subdict.GetNumberOfKeys(), 4);
    QuESo_CHECK_EQUAL( KeySetInfoTypeSubDict::StaticGetAllKeyNames(), "['zero', 'one', 'two', 'three']" );
    BOOST_REQUIRE_THROW( key_set_info_subdict.pGetKey("twelve"), std::exception );

    /// TestKeys4 :: ListTypeTag
    CheckKey<queso::key::ListTypeTag>(TestKeys4::five, 0, "five");
    CheckKey<queso::key::ListTypeTag>(TestKeys4::six, 1, "six");
    CheckKey<queso::key::ListTypeTag>(TestKeys4::seven, 2, "seven");

    typedef decltype(TestKeys4::five)::KeySetInfoType KeySetInfoTypeList;
    KeySetInfoTypeList key_set_info_list{};
    QuESo_CHECK_EQUAL( key_set_info_list.GetNumberOfKeys(), 3);
    QuESo_CHECK_EQUAL( KeySetInfoTypeList::StaticGetAllKeyNames(), "['five', 'six', 'seven']" );
    BOOST_REQUIRE_THROW( key_set_info_list.pGetKey("twelve"), std::exception );
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys5) {
    QuESo_INFO << "Testing :: Test Keys :: Register Triple Key Set" << std::endl;

    /// TestKeys5 :: SubDictTypeTag
    CheckKey<queso::key::SubDictTypeTag>(TestKeys5::zero, 0, "zero");
    CheckKey<queso::key::SubDictTypeTag>(TestKeys5::one, 1, "one");

    typedef decltype(TestKeys5::zero)::KeySetInfoType KeySetInfoTypeSubDict;
    KeySetInfoTypeSubDict key_set_info_subdict{};
    QuESo_CHECK_EQUAL( key_set_info_subdict.GetNumberOfKeys(), 2);
    QuESo_CHECK_EQUAL( KeySetInfoTypeSubDict::StaticGetAllKeyNames(), "['zero', 'one']" );
    BOOST_REQUIRE_THROW( key_set_info_subdict.pGetKey("twelve"), std::exception );

    /// TestKeys5 :: ListTypeTag
    CheckKey<queso::key::ListTypeTag>(TestKeys5::five, 0, "five");
    CheckKey<queso::key::ListTypeTag>(TestKeys5::six, 1, "six");

    typedef decltype(TestKeys5::five)::KeySetInfoType KeySetInfoTypeList;
    KeySetInfoTypeList key_set_info_list{};
    QuESo_CHECK_EQUAL( key_set_info_list.GetNumberOfKeys(), 2);
    QuESo_CHECK_EQUAL( KeySetInfoTypeList::StaticGetAllKeyNames(), "['five', 'six']" );
    BOOST_REQUIRE_THROW( key_set_info_list.pGetKey("twelve"), std::exception );

    /// TestKeys5 :: ListTypeTag
    CheckKey<queso::key::MainValuesTypeTag, double>(TestKeys5::seven, 0, "seven");
    CheckKey<queso::key::MainValuesTypeTag, IndexType>(TestKeys5::eight, 1, "eight");
    CheckKey<queso::key::MainValuesTypeTag, PointType>(TestKeys5::nine, 2, "nine");
    CheckKey<queso::key::MainValuesTypeTag, IndexType>(TestKeys5::ten, 3, "ten");

    typedef decltype(TestKeys5::nine)::KeySetInfoType KeySetInfoTypeValue;
    KeySetInfoTypeValue key_set_info_value{};
    QuESo_CHECK_EQUAL( key_set_info_value.GetNumberOfKeys(), 4);
    QuESo_CHECK_EQUAL( KeySetInfoTypeValue::StaticGetAllKeyNames(), "['seven', 'eight', 'nine', 'ten']" );
    BOOST_REQUIRE_THROW( key_set_info_value.pGetKey("one"), std::exception );
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso