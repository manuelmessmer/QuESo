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
#include <doctest/doctest.h>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/containers/dictionary.hpp"

namespace queso::testing {

    TEST_SUITE("DictionaryTestSuite")
    {

        TEST_CASE("DictionaryTestSuite - Dictionary Key Access Type Test")
        {


            // Create dictionary
            using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
            DictionaryType test_dict(
                DictionaryType::KeySetInfosTypeTag<
                    key::detail::TestKeys5MainValuesTypeTagKeySetInfo,
                    key::detail::TestKeys5SubDictTypeTagKeySetInfo,
                    key::detail::TestKeys5ListTypeTagKeySetInfo>{}
            );

            // Add some values
            CHECK(!test_dict.IsSet(TestKeys5::seven));
            CHECK(!test_dict.IsSet(TestKeys5::eight));
            CHECK(!test_dict.IsSet(TestKeys5::nine));
            CHECK(!test_dict.IsSet(TestKeys5::ten));

            test_dict.SetValue(TestKeys5::seven, 2.0);
            test_dict.SetValue(TestKeys5::eight, 3u);

            /*---- Fill dictionary ----*/
            {
                // Add SubDictionary 1.
                auto p_subdict_1 = MakeUnique<DictionaryType>(DictionaryType::KeySetInfosTypeTag<
                                                              DictionaryType::EmptyKeySetType,
                                                              key::detail::TestKeys1SubDictTypeTagKeySetInfo,
                                                              DictionaryType::EmptyKeySetType>{});

                if constexpr (!NOTDEBUG) {
                    // Dictionary has an empty data set.
                    REQUIRE_THROWS_AS(p_subdict_1->IsSet(TestKeys3::zero), std::exception);
                    REQUIRE_THROWS_AS(
                        p_subdict_1->SetValue(TestKeys3::zero, PointType({ 2.0, 3.2, 1.1 })), std::exception
                    );
                }

                {  // Add SubDictionary 2.
                    auto p_subdict_2 = MakeUnique<DictionaryType>(DictionaryType::KeySetInfosTypeTag<
                                                                  key::detail::TestKeys3MainValuesTypeTagKeySetInfo,
                                                                  DictionaryType::EmptyKeySetType,
                                                                  DictionaryType::EmptyKeySetType>{});

                    if constexpr (!NOTDEBUG) {
                        // Dictionary does not contain any subdictionaries.
                        REQUIRE_THROWS_AS(
                            p_subdict_2->SetSubDictionary(TestKeys1::zero, std::move(p_subdict_2)), std::exception
                        );
                    }

                    p_subdict_2->SetValue(TestKeys3::zero, PointType({ 2.0, 3.2, 1.1 }));
                    p_subdict_2->SetValue(TestKeys3::one, Vector3i({ 1, 2, 3 }));
                    p_subdict_2->SetValue(TestKeys3::two, false);
                    p_subdict_2->SetValue(TestKeys3::four, 4u);
                    p_subdict_2->SetValue(TestKeys3::seven, GridType::b_spline_grid);

                    if constexpr (!NOTDEBUG) {
                        // Wrong key type.
                        REQUIRE_THROWS_AS(
                            p_subdict_1->SetSubDictionary(TestKeys4::zero, std::move(p_subdict_2)), std::exception
                        );
                    }
                    p_subdict_1->SetSubDictionary(TestKeys1::zero, std::move(p_subdict_2));
                }

                test_dict.SetSubDictionary(TestKeys5::zero, std::move(p_subdict_1));

                // Add List  containig SubDictionary 3.
                auto p_subdict_3 = MakeUnique<DictionaryType>(DictionaryType::KeySetInfosTypeTag<
                                                              key::detail::TestKeys3MainValuesTypeTagKeySetInfo,
                                                              DictionaryType::EmptyKeySetType,
                                                              DictionaryType::EmptyKeySetType>{});

                p_subdict_3->SetValue(TestKeys3::zero, PointType({ 4.0, 5.0, 6.0 }));
                p_subdict_3->SetValue(TestKeys3::three, 3.0);
                p_subdict_3->SetValue(TestKeys3::four, 3u);
                p_subdict_3->SetValue(TestKeys3::five, std::string("Hello"));
                p_subdict_3->SetValue(TestKeys3::seven, GridType::hexahedral_fe_grid);

                if constexpr (!NOTDEBUG) {
                    // Dictionary does not contain any list.
                    REQUIRE_THROWS_AS(p_subdict_3->GetList(TestKeys5::five), std::exception);
                    // Wrong key type.
                    REQUIRE_THROWS_AS(test_dict.GetList(TestKeys4::five), std::exception);
                }
                auto& list_a = test_dict.GetList(TestKeys5::five);
                list_a.push_back(std::move(p_subdict_3));
            }

            /*---- Check dictionary ----*/

            // Values are set.
            CHECK(test_dict.IsSet(TestKeys5::seven));
            CHECK(test_dict.IsSet(TestKeys5::eight));

            CheckNear(test_dict.TryGetValue<double>(TestKeys5::seven)->get(), 2.0, 1e-10);
            CHECK((test_dict.TryGetValue<IndexType>(TestKeys5::eight)->get()) == (3u));

            // Values are not set.
            CHECK(!test_dict.IsSet(TestKeys5::nine));
            CHECK(!test_dict.IsSet(TestKeys5::ten));

            // Check Subdictionary.
            auto& r_subdict_1 = test_dict[TestKeys5::zero];

            if constexpr (!NOTDEBUG) {
                // Dictionary has an empty data set.
                REQUIRE_THROWS_AS(r_subdict_1.IsSet(TestKeys3::zero), std::exception);
                REQUIRE_THROWS_AS(r_subdict_1.TryGetValue<Vector3i>(TestKeys3::one), std::exception);
                // Wrong key.
                REQUIRE_THROWS_AS(r_subdict_1[TestKeys4::zero], std::exception);
                // Dict not set.
                REQUIRE_THROWS_AS(r_subdict_1[TestKeys1::one], std::exception);
            }

            const auto& r_subdict_2 = r_subdict_1[TestKeys1::zero];
            if constexpr (!NOTDEBUG) {
                // Dictionary does not contain any subdictionaries.
                REQUIRE_THROWS_AS(r_subdict_2[TestKeys1::zero], std::exception);
            }

            CHECK(r_subdict_2.IsSet(TestKeys3::zero));
            CHECK(r_subdict_2.IsSet(TestKeys3::one));
            CHECK(r_subdict_2.IsSet(TestKeys3::two));
            CHECK(r_subdict_2.IsSet(TestKeys3::four));
            CHECK(r_subdict_2.IsSet(TestKeys3::seven));
            r_subdict_2.CheckRequired();

            CheckPointNear(
                r_subdict_2.GetRequiredValue<PointType>(TestKeys3::zero), PointType({ 2.0, 3.2, 1.1 }), 1e-10
            );
            CHECK(r_subdict_2.TryGetValue<Vector3i>(TestKeys3::one)->get() == Vector3i({ 1, 2, 3 }));
            CHECK(!r_subdict_2.TryGetValue<bool>(TestKeys3::two)->get());
            CHECK((r_subdict_2.GetRequiredValue<IndexType>(TestKeys3::four)) == (4u));
            CHECK((r_subdict_2.GetRequiredValue<GridTypeType>(TestKeys3::seven)) == (GridType::b_spline_grid));

            CHECK(!r_subdict_2.IsSet(TestKeys3::five));

            // Check List (Note: list has just one item).
            for (const auto& items : test_dict.GetList(TestKeys5::five)) {
                items->CheckRequired();
                CheckPointNear(
                    items->GetRequiredValue<PointType>(TestKeys3::zero), PointType({ 4.0, 5.0, 6.0 }), 1e-10
                );
                CheckNear(items->TryGetValue<double>(TestKeys3::three)->get(), 3.0, 1e-10);
                CHECK((items->GetRequiredValue<IndexType>(TestKeys3::four)) == (3u));
                CHECK((items->TryGetValue<std::string>(TestKeys3::five)->get()) == ("Hello"));
                CHECK((items->GetRequiredValue<GridTypeType>(TestKeys3::seven)) == (GridType::hexahedral_fe_grid));
            }
        }

        TEST_CASE("DictionaryTestSuite - Dictionary String Access Type Test")
        {


            // Create dictionary
            using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
            DictionaryType test_dict(
                DictionaryType::KeySetInfosTypeTag<
                    key::detail::TestKeys5MainValuesTypeTagKeySetInfo,
                    key::detail::TestKeys5SubDictTypeTagKeySetInfo,
                    key::detail::TestKeys5ListTypeTagKeySetInfo>{}
            );

            using StringAccess = DictionaryStringAccess<DictionaryType>;

            // Check has
            CHECK(StringAccess::Has(test_dict, "zero"));
            CHECK(StringAccess::Has(test_dict, "one"));
            CHECK(StringAccess::Has(test_dict, "five"));
            CHECK(StringAccess::Has(test_dict, "six"));
            CHECK(StringAccess::Has(test_dict, "seven"));

            CHECK(!StringAccess::Has(test_dict, "two"));

            // Add some values
            CHECK(!StringAccess::IsSet(test_dict, "seven"));
            CHECK(!StringAccess::IsSet(test_dict, "eight"));
            CHECK(!StringAccess::IsSet(test_dict, "nine"));
            CHECK(!StringAccess::IsSet(test_dict, "ten"));

            StringAccess::SetValue(test_dict, "seven", 2.0);
            StringAccess::SetValue(test_dict, "eight", 3u);

            /*---- Fill dictionary ----*/
            {
                // Add SubDictionary 1.
                auto p_subdict_1 = MakeUnique<DictionaryType>(DictionaryType::KeySetInfosTypeTag<
                                                              DictionaryType::EmptyKeySetType,
                                                              key::detail::TestKeys1SubDictTypeTagKeySetInfo,
                                                              DictionaryType::EmptyKeySetType>{});

                // Dictionary has an empty data set.
                REQUIRE_THROWS_AS(StringAccess::IsSet(*p_subdict_1, "zero"), std::exception);
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(*p_subdict_1, "zero", PointType({ 2.0, 3.2, 1.1 })), std::exception
                );

                {  // Add SubDictionary 2.
                    auto p_subdict_2 = MakeUnique<DictionaryType>(DictionaryType::KeySetInfosTypeTag<
                                                                  key::detail::TestKeys3MainValuesTypeTagKeySetInfo,
                                                                  DictionaryType::EmptyKeySetType,
                                                                  DictionaryType::EmptyKeySetType>{});

                    // Dictionary does not contain any subdictionaries.
                    REQUIRE_THROWS_AS(
                        StringAccess::SetSubDictionary(*p_subdict_2, "zero", std::move(p_subdict_2)), std::exception
                    );

                    StringAccess::SetValue(*p_subdict_2, "zero", PointType({ 2.0, 3.2, 1.1 }));
                    StringAccess::SetValue(*p_subdict_2, "one", Vector3i({ 1, 2, 3 }));
                    StringAccess::SetValue(*p_subdict_2, "two", false);

                    // Wrong key type.
                    REQUIRE_THROWS_AS(
                        StringAccess::SetSubDictionary(*p_subdict_1, "wrong_key", std::move(p_subdict_2)),
                        std::exception
                    );

                    StringAccess::SetSubDictionary(*p_subdict_1, "zero", std::move(p_subdict_2));
                }

                StringAccess::SetSubDictionary(test_dict, "zero", std::move(p_subdict_1));

                // Add List  containig SubDictionary 3.
                auto p_subdict_3 = MakeUnique<DictionaryType>(DictionaryType::KeySetInfosTypeTag<
                                                              key::detail::TestKeys3MainValuesTypeTagKeySetInfo,
                                                              DictionaryType::EmptyKeySetType,
                                                              DictionaryType::EmptyKeySetType>{});

                StringAccess::SetValue(*p_subdict_3, "three", 3.0);
                StringAccess::SetValue(*p_subdict_3, "four", 3u);
                StringAccess::SetValue(*p_subdict_3, "five", std::string("Hello"));


                // Dictionary does not contain any list.
                REQUIRE_THROWS_AS(StringAccess::GetList(*p_subdict_3, "five"), std::exception);
                // Wrong key type.
                REQUIRE_THROWS_AS(StringAccess::GetList(test_dict, "wrong_key"), std::exception);

                auto& list_a = StringAccess::GetList(test_dict, "five");
                list_a.push_back(std::move(p_subdict_3));
            }

            /*---- Check dictionary ----*/

            // Values are set.
            CHECK(StringAccess::IsSet(test_dict, "seven"));
            CHECK(StringAccess::IsSet(test_dict, "eight"));

            CheckNear(StringAccess::GetValue<double>(test_dict, "seven"), 2.0, 1e-10);
            CHECK((StringAccess::GetValue<IndexType>(test_dict, "eight")) == (3u));

            // Values are not set.
            CHECK(!StringAccess::IsSet(test_dict, "nine"));
            CHECK(!StringAccess::IsSet(test_dict, "ten"));

            // Check Subdictionary.
            auto& r_subdict_1 = StringAccess::GetSubDictionary(test_dict, "zero");

            // Dictionary has an empty data set.
            REQUIRE_THROWS_AS(StringAccess::IsSet(r_subdict_1, "zero"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<PointType>(r_subdict_1, "zero"), std::exception);
            // Wrong key.
            REQUIRE_THROWS_AS(StringAccess::GetSubDictionary(r_subdict_1, "wrong_key"), std::exception);
            // Dict not set.
            REQUIRE_THROWS_AS(StringAccess::GetSubDictionary(r_subdict_1, "one"), std::exception);

            auto& r_subdict_2 = StringAccess::GetSubDictionary(r_subdict_1, "zero");

            // Dictionary does not contain any subdictionaries.
            REQUIRE_THROWS_AS(StringAccess::GetSubDictionary(r_subdict_2, "zero"), std::exception);

            CHECK(StringAccess::IsSet(r_subdict_2, "zero"));
            CHECK(StringAccess::IsSet(r_subdict_2, "one"));
            CHECK(StringAccess::IsSet(r_subdict_2, "two"));

            CheckPointNear(StringAccess::GetValue<PointType>(r_subdict_2, "zero"), PointType({ 2.0, 3.2, 1.1 }), 1e-10);
            CHECK(StringAccess::GetValue<Vector3i>(r_subdict_2, "one") == Vector3i({ 1, 2, 3 }));
            CHECK(!StringAccess::GetValue<bool>(r_subdict_2, "two"));

            CHECK(!StringAccess::IsSet(r_subdict_2, "five"));

            // Check List (Note: list has just one item).
            for (const auto& items : test_dict.GetList(TestKeys5::five)) {
                CheckNear(StringAccess::GetValue<double>(*items, "three"), 3.0, 1e-10);
                CHECK((StringAccess::GetValue<IndexType>(*items, "four")) == (3u));
                CHECK((StringAccess::GetValue<std::string>(*items, "five")) == ("Hello"));
            }
        }

    }  // TEST_SUITE

} // namespace queso::testing
