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
#include "queso/containers/dictionary.hpp"

namespace queso {

namespace Testing {

BOOST_AUTO_TEST_SUITE( DictionaryTestSuite )

BOOST_AUTO_TEST_CASE(DictionaryKeyAccessTypeTest) {
    QuESo_INFO << "Testing :: Test Dictionary :: Key access" << std::endl;

    // Create dictionary
    using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
    DictionaryType test_dict( DictionaryType::KeySetInfosTypeTag<
        key::detail::TestKeys5MainValuesTypeTagKeySetInfo,
        key::detail::TestKeys5SubDictTypeTagKeySetInfo,
        key::detail::TestKeys5ListTypeTagKeySetInfo>{} );

    // Add some values
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::seven) );
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::eight) );
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::nine) );
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::ten) );

    test_dict.SetValue(TestKeys5::seven, 2.0);
    test_dict.SetValue(TestKeys5::eight, 3u );

    /*---- Fill dictionary ----*/
    {
        // Add SubDictionary 1.
        auto p_subdict_1 = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
            DictionaryType::EmptyKeySetType,
            key::detail::TestKeys1SubDictTypeTagKeySetInfo,
            DictionaryType::EmptyKeySetType>{} );

        if constexpr (!NOTDEBUG) {
            // Dictionary has an empty data set.
            BOOST_REQUIRE_THROW( p_subdict_1->IsSet(TestKeys3::zero), std::exception );
            BOOST_REQUIRE_THROW( p_subdict_1->SetValue(TestKeys3::zero, PointType({2.0, 3.2, 1.1})), std::exception );
        }

        { // Add SubDictionary 2.
            auto p_subdict_2 = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
                key::detail::TestKeys3MainValuesTypeTagKeySetInfo,
                DictionaryType::EmptyKeySetType,
                DictionaryType::EmptyKeySetType>{} );

            if constexpr (!NOTDEBUG) {
                // Dictionary does not contain any subdictionaries.
                BOOST_REQUIRE_THROW( p_subdict_2->SetSubDictionary(TestKeys1::zero, std::move(p_subdict_2)), std::exception );
            }

            p_subdict_2->SetValue(TestKeys3::zero, PointType({2.0, 3.2, 1.1}));
            p_subdict_2->SetValue(TestKeys3::one, Vector3i({1, 2, 3}));
            p_subdict_2->SetValue(TestKeys3::two, false);

            if constexpr (!NOTDEBUG) {
                // Wrong key type.
                BOOST_REQUIRE_THROW( p_subdict_1->SetSubDictionary(TestKeys4::zero, std::move(p_subdict_2)), std::exception );
            }
            p_subdict_1->SetSubDictionary(TestKeys1::zero, std::move(p_subdict_2));
        }

        test_dict.SetSubDictionary(TestKeys5::zero, std::move(p_subdict_1) );

        // Add List  containig SubDictionary 3.
        auto p_subdict_3 = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
            key::detail::TestKeys3MainValuesTypeTagKeySetInfo,
            DictionaryType::EmptyKeySetType,
            DictionaryType::EmptyKeySetType>{} );

        p_subdict_3->SetValue(TestKeys3::three, 3.0);
        p_subdict_3->SetValue(TestKeys3::four, 3u);
        p_subdict_3->SetValue(TestKeys3::five, std::string("Hello"));

        if constexpr (!NOTDEBUG) {
            // Dictionary does not contain any list.
            BOOST_REQUIRE_THROW( p_subdict_3->GetList(TestKeys5::five), std::exception );
            // Wrong key type.
            BOOST_REQUIRE_THROW( test_dict.GetList(TestKeys4::five), std::exception );
        }
        auto& list_a = test_dict.GetList(TestKeys5::five);
        list_a.push_back(std::move(p_subdict_3));
    }

    /*---- Check dictionary ----*/

    // Values are set.
    QuESo_CHECK( test_dict.IsSet(TestKeys5::seven) );
    QuESo_CHECK( test_dict.IsSet(TestKeys5::eight) );

    QuESo_CHECK_NEAR( test_dict.GetValue<double>(TestKeys5::seven), 2.0, 1e-10);
    QuESo_CHECK_EQUAL( test_dict.GetValue<IndexType>(TestKeys5::eight ), 3u );
    QuESo_CHECK_NEAR( test_dict.GetValueFast<double>(TestKeys5::seven), 2.0, 1e-10);
    QuESo_CHECK_EQUAL( test_dict.GetValueFast<IndexType>(TestKeys5::eight ), 3u );

    // Values are not set.
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::nine) );
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::ten) );

    // Check Subdictionary.
    auto& r_subdict_1 = test_dict[TestKeys5::zero];

    if constexpr (!NOTDEBUG) {
        // Dictionary has an empty data set.
        BOOST_REQUIRE_THROW( r_subdict_1.IsSet(TestKeys3::zero), std::exception );
        BOOST_REQUIRE_THROW( r_subdict_1.GetValue<PointType>(TestKeys3::zero), std::exception );
        // Wrong key.
        BOOST_REQUIRE_THROW( r_subdict_1[TestKeys4::zero], std::exception );
        // Dict not set.
        BOOST_REQUIRE_THROW( r_subdict_1[TestKeys1::one], std::exception );
    }

    const auto& r_subdict_2 = r_subdict_1[TestKeys1::zero];
    if constexpr (!NOTDEBUG) {
        // Dictionary does not contain any subdictionaries.
        BOOST_REQUIRE_THROW( r_subdict_2[TestKeys1::zero], std::exception );

    }

    QuESo_CHECK( r_subdict_2.IsSet(TestKeys3::zero) );
    QuESo_CHECK( r_subdict_2.IsSet(TestKeys3::one) );
    QuESo_CHECK( r_subdict_2.IsSet(TestKeys3::two) );

    QuESo_CHECK_POINT_NEAR( r_subdict_2.GetValue<PointType>(TestKeys3::zero), PointType({2.0, 3.2, 1.1}), 1e-10 );
    QuESo_CHECK_Vector3i_EQUAL( r_subdict_2.GetValue<Vector3i>(TestKeys3::one), Vector3i({1, 2, 3}) );
    QuESo_CHECK( !r_subdict_2.GetValue<bool>(TestKeys3::two) );

    QuESo_CHECK( !r_subdict_2.IsSet(TestKeys3::five) );

    // Check List (Note: list has just one item).
    for( const auto& items : test_dict.GetList(TestKeys5::five) ) {
        QuESo_CHECK_NEAR( items->GetValue<double>(TestKeys3::three), 3.0, 1e-10);
        QuESo_CHECK_EQUAL( items->GetValue<IndexType>(TestKeys3::four), 3u);
        QuESo_CHECK_EQUAL( items->GetValue<std::string>(TestKeys3::five), "Hello");

        QuESo_CHECK_NEAR( items->GetValueFast<double>(TestKeys3::three), 3.0, 1e-10);
        QuESo_CHECK_EQUAL( items->GetValueFast<IndexType>(TestKeys3::four), 3u);
        QuESo_CHECK_EQUAL( items->GetValueFast<std::string>(TestKeys3::five), "Hello");
    }
}

BOOST_AUTO_TEST_CASE(DictionaryStringAccessTypeTest) {
    QuESo_INFO << "Testing :: Test Dictionary :: String access" << std::endl;

    // Create dictionary
    using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
    DictionaryType test_dict( DictionaryType::KeySetInfosTypeTag<
        key::detail::TestKeys5MainValuesTypeTagKeySetInfo,
        key::detail::TestKeys5SubDictTypeTagKeySetInfo,
        key::detail::TestKeys5ListTypeTagKeySetInfo>{} );

    using StringAccess = DictionaryStringAccess<DictionaryType>;

    // Check has
    QuESo_CHECK( StringAccess::Has(test_dict, "zero") )
    QuESo_CHECK( StringAccess::Has(test_dict, "one") )
    QuESo_CHECK( StringAccess::Has(test_dict, "five") )
    QuESo_CHECK( StringAccess::Has(test_dict, "six") )
    QuESo_CHECK( StringAccess::Has(test_dict, "seven") )

    QuESo_CHECK( !StringAccess::Has(test_dict, "two") )

    // Add some values
    QuESo_CHECK( !StringAccess::IsSet(test_dict, "seven") );
    QuESo_CHECK( !StringAccess::IsSet(test_dict, "eight") );
    QuESo_CHECK( !StringAccess::IsSet(test_dict, "nine") );
    QuESo_CHECK( !StringAccess::IsSet(test_dict, "ten") );

    StringAccess::SetValue(test_dict, "seven", 2.0);
    StringAccess::SetValue(test_dict, "eight", 3u );

    /*---- Fill dictionary ----*/
    {
        // Add SubDictionary 1.
        auto p_subdict_1 = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
            DictionaryType::EmptyKeySetType,
            key::detail::TestKeys1SubDictTypeTagKeySetInfo,
            DictionaryType::EmptyKeySetType>{} );

        // Dictionary has an empty data set.
        BOOST_REQUIRE_THROW( StringAccess::IsSet(*p_subdict_1, "zero"), std::exception );
        BOOST_REQUIRE_THROW( StringAccess::SetValue(*p_subdict_1, "zero", PointType({2.0, 3.2, 1.1})), std::exception );

        { // Add SubDictionary 2.
            auto p_subdict_2 = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
                key::detail::TestKeys3MainValuesTypeTagKeySetInfo,
                DictionaryType::EmptyKeySetType,
                DictionaryType::EmptyKeySetType>{} );

            // Dictionary does not contain any subdictionaries.
            BOOST_REQUIRE_THROW( StringAccess::SetSubDictionary(*p_subdict_2, "zero", std::move(p_subdict_2)), std::exception );

            StringAccess::SetValue(*p_subdict_2, "zero", PointType({2.0, 3.2, 1.1}));
            StringAccess::SetValue(*p_subdict_2, "one", Vector3i({1, 2, 3}));
            StringAccess::SetValue(*p_subdict_2, "two", false);

            // Wrong key type.
            BOOST_REQUIRE_THROW( StringAccess::SetSubDictionary(*p_subdict_1, "wrong_key", std::move(p_subdict_2)), std::exception );

            StringAccess::SetSubDictionary(*p_subdict_1, "zero", std::move(p_subdict_2));
        }

        StringAccess::SetSubDictionary(test_dict, "zero", std::move(p_subdict_1) );

        // Add List  containig SubDictionary 3.
        auto p_subdict_3 = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
            key::detail::TestKeys3MainValuesTypeTagKeySetInfo,
            DictionaryType::EmptyKeySetType,
            DictionaryType::EmptyKeySetType>{} );

        StringAccess::SetValue(*p_subdict_3, "three", 3.0);
        StringAccess::SetValue(*p_subdict_3, "four", 3u);
        StringAccess::SetValue(*p_subdict_3, "five", std::string("Hello"));


        // Dictionary does not contain any list.
        BOOST_REQUIRE_THROW( StringAccess::GetList(*p_subdict_3, "five"), std::exception );
        // Wrong key type.
        BOOST_REQUIRE_THROW( StringAccess::GetList(test_dict, "wrong_key"), std::exception );

        auto& list_a = StringAccess::GetList(test_dict, "five");
        list_a.push_back(std::move(p_subdict_3));
    }

    /*---- Check dictionary ----*/

    // Values are set.
    QuESo_CHECK( StringAccess::IsSet(test_dict, "seven") );
    QuESo_CHECK( StringAccess::IsSet(test_dict, "eight") );

    QuESo_CHECK_NEAR( StringAccess::GetValue<double>(test_dict, "seven"), 2.0, 1e-10);
    QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(test_dict, "eight" ), 3u );

    // Values are not set.
    QuESo_CHECK( !StringAccess::IsSet(test_dict, "nine") );
    QuESo_CHECK( !StringAccess::IsSet(test_dict, "ten") );

    // Check Subdictionary.
    auto& r_subdict_1 = StringAccess::GetSubDictionary(test_dict, "zero");

    // Dictionary has an empty data set.
    BOOST_REQUIRE_THROW( StringAccess::IsSet(r_subdict_1, "zero"), std::exception );
    BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(r_subdict_1, "zero"), std::exception );
    // Wrong key.
    BOOST_REQUIRE_THROW( StringAccess::GetSubDictionary(r_subdict_1, "wrong_key"), std::exception );
    // Dict not set.
    BOOST_REQUIRE_THROW( StringAccess::GetSubDictionary(r_subdict_1, "one"), std::exception );

    auto& r_subdict_2 = StringAccess::GetSubDictionary(r_subdict_1, "zero");

    // Dictionary does not contain any subdictionaries.
    BOOST_REQUIRE_THROW( StringAccess::GetSubDictionary(r_subdict_2, "zero"), std::exception );

    QuESo_CHECK( StringAccess::IsSet(r_subdict_2, "zero") );
    QuESo_CHECK( StringAccess::IsSet(r_subdict_2, "one") );
    QuESo_CHECK( StringAccess::IsSet(r_subdict_2, "two") );

    QuESo_CHECK_POINT_NEAR( StringAccess::GetValue<PointType>(r_subdict_2, "zero"), PointType({2.0, 3.2, 1.1}), 1e-10 );
    QuESo_CHECK_Vector3i_EQUAL( StringAccess::GetValue<Vector3i>(r_subdict_2, "one"), Vector3i({1, 2, 3}) );
    QuESo_CHECK( !StringAccess::GetValue<bool>(r_subdict_2, "two") );

    QuESo_CHECK( !StringAccess::IsSet(r_subdict_2, "five") );

    // Check List (Note: list has just one item).
    for( const auto& items : test_dict.GetList(TestKeys5::five) ) {
        QuESo_CHECK_NEAR( StringAccess::GetValue<double>(*items, "three"), 3.0, 1e-10);
        QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(*items, "four"), 3u);
        QuESo_CHECK_EQUAL( StringAccess::GetValue<std::string>(*items, "five"), "Hello");
    }
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

