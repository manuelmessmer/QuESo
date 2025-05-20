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
#include "queso/containers/dictionary_tmp.hpp"

namespace queso {


namespace Testing {

BOOST_AUTO_TEST_SUITE( DictionaryTestSuite )

BOOST_AUTO_TEST_CASE(SettingsWrongTypeTest) {
    QuESo_INFO << "Testing :: Test Dictionary :: Key access" << std::endl;

    using SettingsBaseType = Dictionary<queso::key::MainValuesTypeTag>;
    SettingsBaseType test_dict( SettingsBaseType::KeySetInfoTypeTag<
        key::detail::TestKeys5MainValuesTypeTagKeySetInfo,
        key::detail::TestKeys5SubDictTypeTagKeySetInfo,
        key::detail::TestKeys5ListTypeTagKeySetInfo>{} );

    // Values
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::seven) );
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::eight) );
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::nine) );
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::ten) );

    test_dict.SetValue(TestKeys5::seven, 2.0);
    test_dict.SetValue(TestKeys5::eight, 3u );

    { // Fill Dictionary
        // SubDictionary 1
        auto p_subdict_1 = MakeUnique<SettingsBaseType>( SettingsBaseType::KeySetInfoTypeTag<
            std::false_type,
            key::detail::TestKeys1SubDictTypeTagKeySetInfo,
            std::false_type>{} );

        if constexpr (!NOTDEBUG) {
            // Dictionary has an empty data set.
            BOOST_REQUIRE_THROW( p_subdict_1->IsSet(TestKeys3::zero), std::exception );
            BOOST_REQUIRE_THROW( p_subdict_1->SetValue(TestKeys3::zero, PointType({2.0, 3.2, 1.1})), std::exception );
        }

        {
            auto p_subdict_2 = MakeUnique<SettingsBaseType>( SettingsBaseType::KeySetInfoTypeTag<
                key::detail::TestKeys3MainValuesTypeTagKeySetInfo,
                std::false_type,
                std::false_type>{} );

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

        // Add a dict to a list
        auto p_subdict_3 = MakeUnique<SettingsBaseType>( SettingsBaseType::KeySetInfoTypeTag<
            key::detail::TestKeys3MainValuesTypeTagKeySetInfo,
            std::false_type,
            std::false_type>{} );

        p_subdict_3->SetValue(TestKeys3::three, 3.0);
        p_subdict_3->SetValue(TestKeys3::four, 3u);
        p_subdict_3->SetValue(TestKeys3::five, std::string("Hello"));

        if constexpr (!NOTDEBUG) {
            // Dictionary does not contain any list.
            BOOST_REQUIRE_THROW( p_subdict_3->GetList(TestKeys5::five), std::exception );
            // Wrong key type
            BOOST_REQUIRE_THROW( test_dict.GetList(TestKeys4::five), std::exception );
        }
        auto& list_a = test_dict.GetList(TestKeys5::five);
        list_a.push_back(std::move(p_subdict_3));
    }

    /* Check Dictionary */

    // Values are set
    QuESo_CHECK( test_dict.IsSet(TestKeys5::seven) );
    QuESo_CHECK( test_dict.IsSet(TestKeys5::eight) );

    QuESo_CHECK_NEAR( test_dict.GetValue<double>(TestKeys5::seven), 2.0, 1e-10);
    QuESo_CHECK_EQUAL( test_dict.GetValue<IndexType>(TestKeys5::eight ), 3u );
    QuESo_CHECK_NEAR( test_dict.GetValueFast<double>(TestKeys5::seven), 2.0, 1e-10);
    QuESo_CHECK_EQUAL( test_dict.GetValueFast<IndexType>(TestKeys5::eight ), 3u );

    // Values are not set
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::nine) );
    QuESo_CHECK( !test_dict.IsSet(TestKeys5::ten) );

    // Check Subdictionary
    auto& r_subdict_1 = test_dict[TestKeys5::zero];

    if constexpr (!NOTDEBUG) {
        // Dictionary has an empty data set.
        BOOST_REQUIRE_THROW( r_subdict_1.IsSet(TestKeys3::zero), std::exception );
        BOOST_REQUIRE_THROW( r_subdict_1.GetValue<PointType>(TestKeys3::zero), std::exception );
        // Wrong key
        BOOST_REQUIRE_THROW( r_subdict_1[TestKeys4::zero], std::exception );
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

    // Check List (Note: list has just one item)
    for( const auto& items : test_dict.GetList(TestKeys5::five) ) {
        QuESo_CHECK_NEAR( items->GetValue<double>(TestKeys3::three), 3.0, 1e-10);
        QuESo_CHECK_EQUAL( items->GetValue<IndexType>(TestKeys3::four), 3u);
        QuESo_CHECK_EQUAL( items->GetValue<std::string>(TestKeys3::five), "Hello");

        QuESo_CHECK_NEAR( items->GetValueFast<double>(TestKeys3::three), 3.0, 1e-10);
        QuESo_CHECK_EQUAL( items->GetValueFast<IndexType>(TestKeys3::four), 3u);
        QuESo_CHECK_EQUAL( items->GetValueFast<std::string>(TestKeys3::five), "Hello");
    }

}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

