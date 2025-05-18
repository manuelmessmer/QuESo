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

template<typename TKeyType, typename TKeySetToWhat, typename TKeyToWhat = TKeySetToWhat>
constexpr bool CheckKeyStatic(const TKeyType& rKey, IndexType Index, const char* rName) {

    QuESo_CHECK_EQUAL( rKey.Index(), Index );
    QuESo_CHECK_EQUAL( rKey.Name(), rName );

    static_assert( std::is_same_v<typename TKeyType::KeyToWhat, TKeyToWhat> );
    static_assert( std::is_same_v<typename TKeyType::KeySetInfoType::KeySetToWhat, TKeySetToWhat> );

    return true;
}

template<typename TKeyType, typename TKeyToWhat>
void CheckKeyDynamic(const TKeyType& rKey, IndexType Index, const char* rName) {

    typename TKeyType::KeySetInfoType key_set_info;
    const auto p_key = key_set_info.pGetKey(rKey.Name());

    QuESo_CHECK(p_key->TargetTypeIndex() == std::type_index(typeid(TKeyToWhat)));
    QuESo_CHECK(p_key->KeySetInfoTypeIndex() == std::type_index(typeid(typename TKeyType::KeySetInfoType)));

    QuESo_CHECK( key_set_info.IsCorrectKeyType(*p_key) );

    QuESo_CHECK_EQUAL(p_key->Name(), rKey.Name());
    QuESo_CHECK_EQUAL(p_key->Index(), rKey.Index());
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys1) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 1" << std::endl;

    static_assert(TestKeys1::zero.Index() == 0);
    static_assert(TestKeys1::zero.Name() == "zero");

    using KeyType = decltype(TestKeys1::zero);

    static_assert( CheckKeyStatic<KeyType, queso::key::SubDictTypeTag>(TestKeys1::zero, 0, "zero") );
    static_assert( CheckKeyStatic<KeyType, queso::key::SubDictTypeTag>(TestKeys1::one, 1, "one") );
    static_assert( CheckKeyStatic<KeyType, queso::key::SubDictTypeTag>(TestKeys1::two, 2, "two") );
    static_assert( CheckKeyStatic<KeyType, queso::key::SubDictTypeTag>(TestKeys1::three, 3, "three") );
    static_assert( CheckKeyStatic<KeyType, queso::key::SubDictTypeTag>(TestKeys1::four, 4, "four") );

    CheckKeyDynamic<KeyType, queso::key::SubDictTypeTag>(TestKeys1::zero, 0, "zero");
    CheckKeyDynamic<KeyType, queso::key::SubDictTypeTag>(TestKeys1::one, 1, "one");
    CheckKeyDynamic<KeyType, queso::key::SubDictTypeTag>(TestKeys1::two, 2, "two");
    CheckKeyDynamic<KeyType, queso::key::SubDictTypeTag>(TestKeys1::three, 3, "three");
    CheckKeyDynamic<KeyType, queso::key::SubDictTypeTag>(TestKeys1::four, 4, "four");

    using KeySetInfoType = KeyType::KeySetInfoType;
    KeySetInfoType key_set_info{};
    QuESo_CHECK_EQUAL( key_set_info.GetNumberOfKeys(), 5);
    QuESo_CHECK_EQUAL( KeySetInfoType::StaticGetAllKeyNames(), "['zero', 'one', 'two', 'three', 'four']" );
    BOOST_REQUIRE_THROW( key_set_info.pGetKey("five"), std::exception );
    static_assert( queso::key::SubDictTypeTag::GetTypeTagName() == "SubDictTypeTag" );

    using KeySetInfoTypeWrong = decltype(TestKeys2::zero)::KeySetInfoType;
    KeySetInfoTypeWrong key_set_info_wrong{};
    auto p_wrong_key = key_set_info_wrong.pGetKey("zero");
    QuESo_CHECK( !key_set_info.IsCorrectKeyType(*p_wrong_key) );
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys2) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 2" << std::endl;

    using KeyType = decltype(TestKeys2::zero);

    /// TestKeys2
    static_assert( CheckKeyStatic<KeyType, queso::key::ListTypeTag>(TestKeys2::zero, 0, "zero") );
    static_assert( CheckKeyStatic<KeyType, queso::key::ListTypeTag>(TestKeys2::one, 1, "one") );
    static_assert( CheckKeyStatic<KeyType, queso::key::ListTypeTag>(TestKeys2::two, 2, "two") );
    static_assert( CheckKeyStatic<KeyType, queso::key::ListTypeTag>(TestKeys2::three, 3, "three") );

    CheckKeyDynamic<KeyType, queso::key::ListTypeTag>(TestKeys2::zero, 0, "zero");
    CheckKeyDynamic<KeyType, queso::key::ListTypeTag>(TestKeys2::one, 1, "one");
    CheckKeyDynamic<KeyType, queso::key::ListTypeTag>(TestKeys2::two, 2, "two");
    CheckKeyDynamic<KeyType, queso::key::ListTypeTag>(TestKeys2::three, 3, "three");

    using KeySetInfoType = KeyType::KeySetInfoType;
    KeySetInfoType key_set_info{};
    QuESo_CHECK_EQUAL( key_set_info.GetNumberOfKeys(), 4);
    QuESo_CHECK_EQUAL( KeySetInfoType::StaticGetAllKeyNames(), "['zero', 'one', 'two', 'three']" );
    BOOST_REQUIRE_THROW( key_set_info.pGetKey("four"), std::exception );

    static_assert( queso::key::ListTypeTag::GetTypeTagName() == "ListTypeTag" );
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys3) {
    QuESo_INFO << "Testing :: Test Keys :: Register Single Key Set 3" << std::endl;

    /// TestKeys3
    static_assert( CheckKeyStatic<decltype(TestKeys3::zero), queso::key::MainValuesTypeTag, PointType>(TestKeys3::zero, 0, "zero") );
    static_assert( CheckKeyStatic<decltype(TestKeys3::one), queso::key::MainValuesTypeTag, Vector3i>(TestKeys3::one, 1, "one") );
    static_assert( CheckKeyStatic<decltype(TestKeys3::two), queso::key::MainValuesTypeTag, bool>(TestKeys3::two, 2, "two") );
    static_assert( CheckKeyStatic<decltype(TestKeys3::three), queso::key::MainValuesTypeTag, double>(TestKeys3::three, 3, "three") );
    static_assert( CheckKeyStatic<decltype(TestKeys3::four), queso::key::MainValuesTypeTag, IndexType>(TestKeys3::four, 4, "four") );
    static_assert( CheckKeyStatic<decltype(TestKeys3::five), queso::key::MainValuesTypeTag, std::string>(TestKeys3::five, 5, "five") );
    static_assert( CheckKeyStatic<decltype(TestKeys3::six), queso::key::MainValuesTypeTag, IntegrationMethodType>(TestKeys3::six, 6, "six") );
    static_assert( CheckKeyStatic<decltype(TestKeys3::seven), queso::key::MainValuesTypeTag, GridTypeType>(TestKeys3::seven, 7, "seven") );

    CheckKeyDynamic<decltype(TestKeys3::zero), PointType>(TestKeys3::zero, 0, "zero");
    CheckKeyDynamic<decltype(TestKeys3::one), Vector3i>(TestKeys3::one, 1, "one");
    CheckKeyDynamic<decltype(TestKeys3::two), bool>(TestKeys3::two, 2, "two");
    CheckKeyDynamic<decltype(TestKeys3::three), double>(TestKeys3::three, 3, "three");
    CheckKeyDynamic<decltype(TestKeys3::four), IndexType>(TestKeys3::four, 4, "four");
    CheckKeyDynamic<decltype(TestKeys3::five), std::string>(TestKeys3::five, 5, "five");
    CheckKeyDynamic<decltype(TestKeys3::six), IntegrationMethodType>(TestKeys3::six, 6, "six");
    CheckKeyDynamic<decltype(TestKeys3::seven), GridTypeType>(TestKeys3::seven, 7, "seven");

    using KeySetInfoType = decltype(TestKeys3::zero)::KeySetInfoType;
    using KeySetInfoType_2 = decltype(TestKeys3::one)::KeySetInfoType;
    // Note that TestKeys3::zero and TestKeys3::one are not of the same type.
    static_assert( !std::is_same_v<decltype(TestKeys3::zero), decltype(TestKeys3::one)> );
    // Check if both lead to the same KeySetInfo.
    static_assert( std::is_same_v<KeySetInfoType, KeySetInfoType_2> );

    KeySetInfoType key_set_info{};
    QuESo_CHECK_EQUAL( key_set_info.GetNumberOfKeys(), 8);
    QuESo_CHECK_EQUAL( KeySetInfoType::StaticGetAllKeyNames(), "['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven']" );
    BOOST_REQUIRE_THROW( key_set_info.pGetKey("twelve"), std::exception );

    using KeyType = decltype(TestKeys3::two);

    static_assert( KeyType::KeySetInfoType::KeySetToWhat::GetTypeTagName() == "MainValuesTypeTag" );
    static_assert( KeyType::KeySetInfoType::KeySetToWhat::GetValueTypeName<PointType>() == "PointType" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<Vector3i>() == "Vector3i" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<bool>() == "bool" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<double>() == "double" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<IndexType>() == "IndexType" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<std::string>() == "std::string" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<IntegrationMethodType>() == "IntegrationMethodType" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<GridTypeType>() == "GridTypeType" );
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys4) {
    QuESo_INFO << "Testing :: Test Keys :: Register Double Key Set" << std::endl;

    /// TestKeys4 :: SubDictTypeTag
    static_assert( CheckKeyStatic<decltype(TestKeys4::zero), queso::key::SubDictTypeTag>(TestKeys4::zero, 0, "zero") );
    static_assert( CheckKeyStatic<decltype(TestKeys4::one), queso::key::SubDictTypeTag>(TestKeys4::one, 1, "one") );
    static_assert( CheckKeyStatic<decltype(TestKeys4::two), queso::key::SubDictTypeTag>(TestKeys4::two, 2, "two") );
    static_assert( CheckKeyStatic<decltype(TestKeys4::three), queso::key::SubDictTypeTag>(TestKeys4::three, 3, "three") );

    CheckKeyDynamic<decltype(TestKeys4::zero), queso::key::SubDictTypeTag>(TestKeys4::zero, 0, "zero");
    CheckKeyDynamic<decltype(TestKeys4::one), queso::key::SubDictTypeTag>(TestKeys4::one, 1, "one");
    CheckKeyDynamic<decltype(TestKeys4::two), queso::key::SubDictTypeTag>(TestKeys4::two, 2, "two");
    CheckKeyDynamic<decltype(TestKeys4::three), queso::key::SubDictTypeTag>(TestKeys4::three, 3, "three");

    using KeySetInfoTypeSubDict = decltype(TestKeys4::zero)::KeySetInfoType;
    static_assert( KeySetInfoTypeSubDict::KeySetToWhat::GetTypeTagName() == "SubDictTypeTag" );

    KeySetInfoTypeSubDict key_set_info_subdict{};
    QuESo_CHECK_EQUAL( key_set_info_subdict.GetNumberOfKeys(), 4);
    QuESo_CHECK_EQUAL( KeySetInfoTypeSubDict::StaticGetAllKeyNames(), "['zero', 'one', 'two', 'three']" );
    BOOST_REQUIRE_THROW( key_set_info_subdict.pGetKey("twelve"), std::exception );

    /// TestKeys4 :: ListTypeTag
    static_assert( CheckKeyStatic<decltype(TestKeys4::five), queso::key::ListTypeTag>(TestKeys4::five, 0, "five") );
    static_assert( CheckKeyStatic<decltype(TestKeys4::six), queso::key::ListTypeTag>(TestKeys4::six, 1, "six") );
    static_assert( CheckKeyStatic<decltype(TestKeys4::seven), queso::key::ListTypeTag>(TestKeys4::seven, 2, "seven") );

    CheckKeyDynamic<decltype(TestKeys4::five), queso::key::ListTypeTag>(TestKeys4::five, 0, "five");
    CheckKeyDynamic<decltype(TestKeys4::six), queso::key::ListTypeTag>(TestKeys4::six, 1, "six");
    CheckKeyDynamic<decltype(TestKeys4::seven), queso::key::ListTypeTag>(TestKeys4::seven, 2, "seven");

    using KeySetInfoTypeList = decltype(TestKeys4::five)::KeySetInfoType;
    static_assert( KeySetInfoTypeList::KeySetToWhat::GetTypeTagName() == "ListTypeTag" );

    KeySetInfoTypeList key_set_info_list{};
    QuESo_CHECK_EQUAL( key_set_info_list.GetNumberOfKeys(), 3);
    QuESo_CHECK_EQUAL( KeySetInfoTypeList::StaticGetAllKeyNames(), "['five', 'six', 'seven']" );
    BOOST_REQUIRE_THROW( key_set_info_list.pGetKey("twelve"), std::exception );
}

BOOST_AUTO_TEST_CASE(TestRegisterKeys5) {
    QuESo_INFO << "Testing :: Test Keys :: Register Triple Key Set" << std::endl;

    /// TestKeys5 :: SubDictTypeTag
    static_assert( CheckKeyStatic<decltype(TestKeys5::zero), queso::key::SubDictTypeTag>(TestKeys5::zero, 0, "zero") );
    static_assert( CheckKeyStatic<decltype(TestKeys5::one), queso::key::SubDictTypeTag>(TestKeys5::one, 1, "one") );

    CheckKeyDynamic<decltype(TestKeys5::zero), queso::key::SubDictTypeTag>(TestKeys5::zero, 0, "zero");
    CheckKeyDynamic<decltype(TestKeys5::one), queso::key::SubDictTypeTag>(TestKeys5::one, 1, "one");

    using KeySetInfoTypeSubDict = decltype(TestKeys5::zero)::KeySetInfoType;
    using KeySetInfoTypeSubDict_2 = decltype(TestKeys5::one)::KeySetInfoType;
    static_assert( std::is_same_v<decltype(TestKeys5::zero), decltype(TestKeys5::one)> );
    static_assert( std::is_same_v<KeySetInfoTypeSubDict, KeySetInfoTypeSubDict_2> );

    static_assert( KeySetInfoTypeSubDict::KeySetToWhat::GetTypeTagName() == "SubDictTypeTag" );

    KeySetInfoTypeSubDict key_set_info_subdict{};
    QuESo_CHECK_EQUAL( key_set_info_subdict.GetNumberOfKeys(), 2);
    QuESo_CHECK_EQUAL( KeySetInfoTypeSubDict::StaticGetAllKeyNames(), "['zero', 'one']" );
    BOOST_REQUIRE_THROW( key_set_info_subdict.pGetKey("twelve"), std::exception );

    /// TestKeys5 :: ListTypeTag
    static_assert( CheckKeyStatic<decltype(TestKeys5::five), queso::key::ListTypeTag>(TestKeys5::five, 0, "five") );
    static_assert( CheckKeyStatic<decltype(TestKeys5::six), queso::key::ListTypeTag>(TestKeys5::six, 1, "six") );

    CheckKeyDynamic<decltype(TestKeys5::five), queso::key::ListTypeTag>(TestKeys5::five, 0, "five");
    CheckKeyDynamic<decltype(TestKeys5::six), queso::key::ListTypeTag>(TestKeys5::six, 1, "six");

    using KeySetInfoTypeList = decltype(TestKeys5::five)::KeySetInfoType;
    using KeySetInfoTypeList_2 = decltype(TestKeys5::six)::KeySetInfoType;
    static_assert( std::is_same_v<decltype(TestKeys5::five), decltype(TestKeys5::six)> );
    static_assert( std::is_same_v<KeySetInfoTypeList, KeySetInfoTypeList_2> );

    static_assert( KeySetInfoTypeList::KeySetToWhat::GetTypeTagName() == "ListTypeTag" );

    KeySetInfoTypeList key_set_info_list{};
    QuESo_CHECK_EQUAL( key_set_info_list.GetNumberOfKeys(), 2);
    QuESo_CHECK_EQUAL( KeySetInfoTypeList::StaticGetAllKeyNames(), "['five', 'six']" );
    BOOST_REQUIRE_THROW( key_set_info_list.pGetKey("twelve"), std::exception );

    /// TestKeys5 :: MainValuesTypeTag
    static_assert( CheckKeyStatic<decltype(TestKeys5::seven), queso::key::MainValuesTypeTag, double>(TestKeys5::seven, 0, "seven") );
    static_assert( CheckKeyStatic<decltype(TestKeys5::eight), queso::key::MainValuesTypeTag, IndexType>(TestKeys5::eight, 1, "eight") );
    static_assert( CheckKeyStatic<decltype(TestKeys5::nine), queso::key::MainValuesTypeTag, PointType>(TestKeys5::nine, 2, "nine") );
    static_assert( CheckKeyStatic<decltype(TestKeys5::ten), queso::key::MainValuesTypeTag, IndexType>(TestKeys5::ten, 3, "ten") );

    CheckKeyDynamic<decltype(TestKeys5::seven), double>(TestKeys5::seven, 0, "seven");
    CheckKeyDynamic<decltype(TestKeys5::eight), IndexType>(TestKeys5::eight, 1, "eight");
    CheckKeyDynamic<decltype(TestKeys5::nine), PointType>(TestKeys5::nine, 2, "nine");
    CheckKeyDynamic<decltype(TestKeys5::ten), IndexType>(TestKeys5::ten, 3, "ten");

    using KeySetInfoTypeValue = decltype(TestKeys5::nine)::KeySetInfoType;
    using KeySetInfoTypeValue_2 = decltype(TestKeys5::ten)::KeySetInfoType;
    static_assert( !std::is_same_v<decltype(TestKeys5::nine), decltype(TestKeys5::ten)> );
    static_assert( std::is_same_v<KeySetInfoTypeValue, KeySetInfoTypeValue_2> );

    static_assert( KeySetInfoTypeValue::KeySetToWhat::GetTypeTagName() == "MainValuesTypeTag" );
    static_assert( KeySetInfoTypeValue::KeySetToWhat::GetValueTypeName<PointType>() == "PointType" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<Vector3i>() == "Vector3i" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<bool>() == "bool" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<double>() == "double" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<IndexType>() == "IndexType" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<std::string>() == "std::string" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<IntegrationMethodType>() == "IntegrationMethodType" );
    static_assert( queso::key::MainValuesTypeTag::GetValueTypeName<GridTypeType>() == "GridTypeType" );

    KeySetInfoTypeValue key_set_info_value{};
    QuESo_CHECK_EQUAL( key_set_info_value.GetNumberOfKeys(), 4);
    QuESo_CHECK_EQUAL( KeySetInfoTypeValue::StaticGetAllKeyNames(), "['seven', 'eight', 'nine', 'ten']" );
    BOOST_REQUIRE_THROW( key_set_info_value.pGetKey("one"), std::exception );
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso