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
#include "queso/containers/data_set.hpp"

namespace queso::testing {

    TEST_SUITE("DataSetTestSuite")
    {

        TEST_CASE("DataSetTestSuite - Test Data Set Type Traits")
        {


            using DataSetType = DataSet<queso::key::MainValuesTypeTag>;

            enum PlainIntEnum { a1, a2, a3 };
            enum class ClassIntEnum { a1, a2, a3 };

            // Check is_index for
            static_assert(DataSetType::is_index_v<IndexType>);
            static_assert(DataSetType::is_index_v<unsigned int>);
            static_assert(DataSetType::is_index_v<unsigned short>);
            static_assert(DataSetType::is_index_v<unsigned long>);
            static_assert(DataSetType::is_index_v<int>);
            static_assert(DataSetType::is_index_v<short>);
            static_assert(DataSetType::is_index_v<long>);

            static_assert(!DataSetType::is_index_v<PlainIntEnum>);
            static_assert(!DataSetType::is_index_v<ClassIntEnum>);
            static_assert(!DataSetType::is_index_v<float>);
            static_assert(!DataSetType::is_index_v<double>);
            static_assert(!DataSetType::is_index_v<bool>);
            static_assert(!DataSetType::is_index_v<char>);
            static_assert(!DataSetType::is_index_v<signed char>);
            static_assert(!DataSetType::is_index_v<unsigned char>);
            static_assert(!DataSetType::is_index_v<std::string>);

            // Check is_unsigned_index
            static_assert(DataSetType::is_unsigned_index_v<IndexType>);
            static_assert(DataSetType::is_unsigned_index_v<unsigned int>);
            static_assert(DataSetType::is_unsigned_index_v<unsigned short>);
            static_assert(DataSetType::is_unsigned_index_v<unsigned long>);

            static_assert(!DataSetType::is_unsigned_index_v<int>);
            static_assert(!DataSetType::is_unsigned_index_v<short>);
            static_assert(!DataSetType::is_unsigned_index_v<long>);
            static_assert(!DataSetType::is_unsigned_index_v<PlainIntEnum>);
            static_assert(!DataSetType::is_unsigned_index_v<ClassIntEnum>);
            static_assert(!DataSetType::is_unsigned_index_v<float>);
            static_assert(!DataSetType::is_unsigned_index_v<double>);
            static_assert(!DataSetType::is_unsigned_index_v<bool>);
            static_assert(!DataSetType::is_unsigned_index_v<char>);
            static_assert(!DataSetType::is_unsigned_index_v<signed char>);
            static_assert(!DataSetType::is_unsigned_index_v<unsigned char>);
            static_assert(!DataSetType::is_unsigned_index_v<std::string>);

            // Check is_scoped_int_enum
            static_assert(DataSetType::is_key_v<decltype(TestKeys3::zero)>);

            static_assert(!DataSetType::is_key_v<ClassIntEnum>);
            static_assert(!DataSetType::is_key_v<PlainIntEnum>);
            static_assert(!DataSetType::is_key_v<int>);
            static_assert(!DataSetType::is_key_v<IndexType>);

            // Runtime dummy check
            CHECK(true);
        }

        TEST_CASE("DataSetTestSuite - Test Data Set Check Required")
        {


            using DataSetType = DataSet<queso::key::MainValuesTypeTag>;

            {
                DataSetType data_set(
                    DataSetType::KeySetInfoTypeTag<queso::key::detail::MainSettingsMainValuesTypeTagKeySetInfo>{}
                );
                REQUIRE_THROWS_AS(data_set.CheckRequired(), std::exception);

                data_set.SetValue(MainSettings::component_name, std::string("main"));
                data_set.SetValue(MainSettings::input_filename, std::string("input.stl"));
                data_set.SetValue(MainSettings::output_directory_name, std::string("output"));
                data_set.SetValue(MainSettings::echo_level, 1u);
                data_set.SetValue(MainSettings::write_output_to_file, true);
                data_set.CheckRequired();
            }

            {
                DataSetType data_set(
                    DataSetType::KeySetInfoTypeTag<
                        queso::key::detail::BackgroundGridSettingsMainValuesTypeTagKeySetInfo>{}
                );
                data_set.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
                data_set.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType({ 0.0, 0.0, 0.0 }));
                data_set.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType({ 1.0, 1.0, 1.0 }));
                data_set.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType({ 0.0, 0.0, 0.0 }));
                data_set.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType({ 1.0, 1.0, 1.0 }));
                data_set.SetValue(BackgroundGridSettings::polynomial_order, Vector3i({ 2, 2, 2 }));
                REQUIRE_THROWS_AS(data_set.CheckRequired(), std::exception);

                data_set.SetValue(BackgroundGridSettings::number_of_elements, Vector3i({ 4, 4, 4 }));
                data_set.CheckRequired();
            }

            {
                DataSetType data_set(
                    DataSetType::KeySetInfoTypeTag<key::detail::TestKeys3MainValuesTypeTagKeySetInfo>{}
                );
                REQUIRE_THROWS_AS(data_set.CheckRequired(), std::exception);

                data_set.SetValue(TestKeys3::zero, PointType{ 1.0, 2.0, 3.0 });
                REQUIRE_THROWS_AS(data_set.CheckRequired(), std::exception);

                data_set.SetValue(TestKeys3::four, IndexType(4));
                REQUIRE_THROWS_AS(data_set.CheckRequired(), std::exception);

                data_set.SetValue(TestKeys3::seven, GridType::b_spline_grid);
                data_set.CheckRequired();
            }
        }

        TEST_CASE("DataSetTestSuite - Test Data Set Value")
        {


            using DataSetType = DataSet<queso::key::MainValuesTypeTag>;
            DataSetType data_set(DataSetType::KeySetInfoTypeTag<key::detail::TestKeys3MainValuesTypeTagKeySetInfo>{});

            data_set.SetValue(TestKeys3::zero, PointType{ 1.0, 2.0, 3.0 });
            data_set.SetValue(TestKeys3::one, Vector3i{ 1, 2, 3 });
            data_set.SetValue(TestKeys3::two, true);
            data_set.SetValue(TestKeys3::three, 2.0);
            data_set.SetValue(TestKeys3::four, (unsigned int)1);
            data_set.SetValue(TestKeys3::four, (unsigned short)1);
            data_set.SetValue(TestKeys3::four, (unsigned long)1);
            data_set.SetValue(TestKeys3::four, IndexType(1));
            data_set.SetValue(TestKeys3::five, std::string("Hallo"));
            data_set.SetValue(TestKeys3::six, IntegrationMethod::gauss);
            data_set.SetValue(TestKeys3::seven, GridType::b_spline_grid);

            if constexpr (!NOTDEBUG) {
                // Wrong keys
                REQUIRE_THROWS_AS(data_set.SetValue(TestKeys5::seven, 2.0), std::exception);
                REQUIRE_THROWS_AS(data_set.SetValue(TestKeys5::eight, IndexType(2)), std::exception);
                REQUIRE_THROWS_AS(data_set.SetValue(TestKeys5::nine, PointType{ 1.0, 2.0, 3.0 }), std::exception);
                REQUIRE_THROWS_AS(data_set.SetValue(TestKeys5::ten, 4ul), std::exception);
            }

            using StringAccess = DataSetStringAccess<DataSetType>;

            StringAccess::SetValue(data_set, "zero", PointType{ 1.0, 2.0, 3.0 });
            StringAccess::SetValue(
                data_set, "zero", Vector3i{ 1, 2, 3 }
            );  // The conversion from Vector3i to PointType is legal.
            REQUIRE_THROWS_AS(StringAccess::SetValue(data_set, "zero", 1.0), std::exception);

            StringAccess::SetValue(data_set, "one", Vector3i{ 1, 2, 3 });
            REQUIRE_THROWS_AS(StringAccess::SetValue(data_set, "one", PointType{ 1.0, 2.0, 3.0 }), std::exception);

            StringAccess::SetValue(data_set, "two", false);
            REQUIRE_THROWS_AS(StringAccess::SetValue(data_set, "two", 1), std::exception);

            StringAccess::SetValue(data_set, "three", 3.0);
            StringAccess::SetValue(data_set, "three", 3);  // The conversion from int to double is legal.
            StringAccess::SetValue(data_set, "three", IndexType(3));  // The conversion from int to double is legal.
            REQUIRE_THROWS_AS(StringAccess::SetValue(data_set, "three", true), std::exception);

            StringAccess::SetValue(data_set, "four", (int)1);
            StringAccess::SetValue(data_set, "four", (unsigned int)1);
            StringAccess::SetValue(data_set, "four", IndexType(1));
            StringAccess::SetValue(data_set, "four", (long)1);
            StringAccess::SetValue(data_set, "four", (short)1);
            REQUIRE_THROWS_AS(StringAccess::SetValue(data_set, "four", -1), std::exception);
            REQUIRE_THROWS_AS(StringAccess::SetValue(data_set, "four", 1.0), std::exception);

            StringAccess::SetValue(data_set, "five", std::string("hallo"));
            REQUIRE_THROWS_AS(StringAccess::SetValue(data_set, "five", 1), std::exception);

            StringAccess::SetValue(data_set, "six", IntegrationMethod::gauss);
            REQUIRE_THROWS_AS(StringAccess::SetValue(data_set, "six", Vector3i{ 0, 0, 0 }), std::exception);

            StringAccess::SetValue(data_set, "seven", GridType::b_spline_grid);
            REQUIRE_THROWS_AS(StringAccess::SetValue(data_set, "seven", 1), std::exception);

            REQUIRE_THROWS_AS(StringAccess::SetValue(data_set, "wrong_key", 1), std::exception);
        }

        TEST_CASE("DataSetTestSuite - Test Data Set Get Value Key")
        {


            using DataSetType = DataSet<queso::key::MainValuesTypeTag>;
            DataSetType data_set(DataSetType::KeySetInfoTypeTag<key::detail::TestKeys3MainValuesTypeTagKeySetInfo>{});

            // Values are not set.
            CHECK(!data_set.IsSet(TestKeys3::zero));
            CHECK(!data_set.IsSet(TestKeys3::one));
            CHECK(!data_set.IsSet(TestKeys3::two));
            CHECK(!data_set.IsSet(TestKeys3::three));
            CHECK(!data_set.IsSet(TestKeys3::four));
            CHECK(!data_set.IsSet(TestKeys3::five));
            CHECK(!data_set.IsSet(TestKeys3::six));
            CHECK(!data_set.IsSet(TestKeys3::seven));

            // Values are not set
            REQUIRE(!data_set.TryGetValue<Vector3i>(TestKeys3::one).has_value());
            REQUIRE(!data_set.TryGetValue<bool>(TestKeys3::two).has_value());
            REQUIRE(!data_set.TryGetValue<double>(TestKeys3::three).has_value());
            REQUIRE(!data_set.TryGetValue<std::string>(TestKeys3::five).has_value());
            REQUIRE(!data_set.TryGetValue<IntegrationMethodType>(TestKeys3::six).has_value());
            if constexpr (!NOTDEBUG) {
                REQUIRE_THROWS_AS(data_set.GetRequiredValue<PointType>(TestKeys3::zero), std::exception);
                REQUIRE_THROWS_AS(data_set.GetRequiredValue<IndexType>(TestKeys3::four), std::exception);
                REQUIRE_THROWS_AS(data_set.GetRequiredValue<GridTypeType>(TestKeys3::seven), std::exception);
            }

            // Set values
            data_set.SetValue(TestKeys3::zero, PointType{ 1.0, 2.0, 3.0 });
            data_set.SetValue(TestKeys3::one, Vector3i{ 1, 2, 3 });
            data_set.SetValue(TestKeys3::two, true);
            data_set.SetValue(TestKeys3::three, 2.0);
            data_set.SetValue(TestKeys3::four, IndexType(5));
            data_set.SetValue(TestKeys3::five, std::string("Hallo"));
            data_set.SetValue(TestKeys3::six, IntegrationMethod::gauss);
            data_set.SetValue(TestKeys3::seven, GridType::b_spline_grid);

            // Values are set.
            CHECK(data_set.IsSet(TestKeys3::zero));
            CHECK(data_set.IsSet(TestKeys3::one));
            CHECK(data_set.IsSet(TestKeys3::two));
            CHECK(data_set.IsSet(TestKeys3::three));
            CHECK(data_set.IsSet(TestKeys3::four));
            CHECK(data_set.IsSet(TestKeys3::five));
            CHECK(data_set.IsSet(TestKeys3::six));
            CHECK(data_set.IsSet(TestKeys3::seven));

            // Wrong keys
            data_set.CheckRequired();
            if constexpr (!NOTDEBUG) {
                REQUIRE_THROWS_AS(data_set.IsSet(TestKeys5::seven), std::exception);

                REQUIRE_THROWS_AS(data_set.TryGetValue<double>(TestKeys5::seven), std::exception);
                REQUIRE_THROWS_AS(data_set.TryGetValue<IndexType>(TestKeys5::eight), std::exception);
                REQUIRE_THROWS_AS(data_set.TryGetValue<PointType>(TestKeys5::nine), std::exception);
                REQUIRE_THROWS_AS(data_set.TryGetValue<IndexType>(TestKeys5::ten), std::exception);
                REQUIRE_THROWS_AS(data_set.GetRequiredValue<IndexType>(TestKeys5::eleven), std::exception);
            }

            // Get values
            CheckPointNear(data_set.GetRequiredValue<PointType>(TestKeys3::zero), PointType({ 1.0, 2.0, 3.0 }), 1e-12);
            CHECK(data_set.TryGetValue<Vector3i>(TestKeys3::one)->get() == Vector3i({ 1, 2, 3 }));
            CHECK((data_set.TryGetValue<bool>(TestKeys3::two)->get()) == (true));
            CheckNear(data_set.TryGetValue<double>(TestKeys3::three)->get(), 2.0, 1e-12);
            CHECK((data_set.GetRequiredValue<IndexType>(TestKeys3::four)) == (5));
            CHECK((data_set.TryGetValue<std::string>(TestKeys3::five)->get()) == (std::string("Hallo")));
            CHECK((data_set.TryGetValue<IntegrationMethodType>(TestKeys3::six)->get()) == (IntegrationMethod::gauss));
            CHECK((data_set.GetRequiredValue<GridTypeType>(TestKeys3::seven)) == (GridType::b_spline_grid));
        }


        TEST_CASE("DataSetTestSuite - Test Data Set Get Value Key Name")
        {


            using DataSetType = DataSet<queso::key::MainValuesTypeTag>;
            DataSetType data_set(DataSetType::KeySetInfoTypeTag<key::detail::TestKeys3MainValuesTypeTagKeySetInfo>{});

            using StringAccess = DataSetStringAccess<DataSetType>;

            // Has keys.
            CHECK(StringAccess::Has(data_set, "zero"));
            CHECK(StringAccess::Has(data_set, "one"));
            CHECK(StringAccess::Has(data_set, "two"));
            CHECK(StringAccess::Has(data_set, "three"));
            CHECK(StringAccess::Has(data_set, "four"));
            CHECK(StringAccess::Has(data_set, "five"));
            CHECK(StringAccess::Has(data_set, "six"));
            CHECK(StringAccess::Has(data_set, "seven"));

            // Does not have key.
            CHECK(!StringAccess::Has(data_set, "wrong_key"));

            // Values are not set.
            CHECK(!StringAccess::IsSet(data_set, "zero"));
            CHECK(!StringAccess::IsSet(data_set, "one"));
            CHECK(!StringAccess::IsSet(data_set, "two"));
            CHECK(!StringAccess::IsSet(data_set, "three"));
            CHECK(!StringAccess::IsSet(data_set, "four"));
            CHECK(!StringAccess::IsSet(data_set, "five"));
            CHECK(!StringAccess::IsSet(data_set, "six"));
            CHECK(!StringAccess::IsSet(data_set, "seven"));

            // Values are not set.
            REQUIRE_THROWS_AS(StringAccess::GetValue<PointType>(data_set, "zero"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<Vector3i>(data_set, "one"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<bool>(data_set, "two"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<double>(data_set, "three"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<IndexType>(data_set, "four"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<std::string>(data_set, "five"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<IntegrationMethodType>(data_set, "six"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<GridTypeType>(data_set, "seven"), std::exception);

            // Set values
            StringAccess::SetValue(data_set, "zero", PointType{ 1.0, 2.0, 3.0 });
            StringAccess::SetValue(data_set, "one", Vector3i{ 1, 2, 3 });
            StringAccess::SetValue(data_set, "two", true);
            StringAccess::SetValue(data_set, "three", 2.0);
            StringAccess::SetValue(data_set, "four", IndexType(5));
            StringAccess::SetValue(data_set, "five", std::string("Hallo"));
            StringAccess::SetValue(data_set, "six", IntegrationMethod::gauss);
            StringAccess::SetValue(data_set, "seven", GridType::b_spline_grid);

            // Values are set.
            CHECK(StringAccess::IsSet(data_set, "zero"));
            CHECK(StringAccess::IsSet(data_set, "one"));
            CHECK(StringAccess::IsSet(data_set, "two"));
            CHECK(StringAccess::IsSet(data_set, "three"));
            CHECK(StringAccess::IsSet(data_set, "four"));
            CHECK(StringAccess::IsSet(data_set, "five"));
            CHECK(StringAccess::IsSet(data_set, "six"));
            CHECK(StringAccess::IsSet(data_set, "seven"));

            // Wrong key
            REQUIRE_THROWS_AS(StringAccess::IsSet(data_set, "zero_"), std::exception);

            // Wrong keys
            REQUIRE_THROWS_AS(StringAccess::GetValue<PointType>(data_set, "zero_"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<Vector3i>(data_set, "one_"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<bool>(data_set, "two_"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<double>(data_set, "three_"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<IndexType>(data_set, "four_"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<std::string>(data_set, "five_"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<IntegrationMethodType>(data_set, "six_"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<GridTypeType>(data_set, "seven_"), std::exception);

            // Wrong variable types
            REQUIRE_THROWS_AS(StringAccess::GetValue<Vector3i>(data_set, "zero"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<bool>(data_set, "one"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<double>(data_set, "two"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<IndexType>(data_set, "three"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<double>(data_set, "four"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<GridTypeType>(data_set, "five"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<std::string>(data_set, "six"), std::exception);
            REQUIRE_THROWS_AS(StringAccess::GetValue<IndexType>(data_set, "seven"), std::exception);

            // Get values
            CheckPointNear(StringAccess::GetValue<PointType>(data_set, "zero"), PointType({ 1.0, 2.0, 3.0 }), 1e-12);
            CHECK(StringAccess::GetValue<Vector3i>(data_set, "one") == Vector3i({ 1, 2, 3 }));
            CHECK((StringAccess::GetValue<bool>(data_set, "two")) == (true));
            CheckNear(StringAccess::GetValue<double>(data_set, "three"), 2.0, 1e-12);
            CHECK((StringAccess::GetValue<IndexType>(data_set, "four")) == (5));
            CHECK((StringAccess::GetValue<std::string>(data_set, "five")) == (std::string("Hallo")));
            CHECK((StringAccess::GetValue<IntegrationMethodType>(data_set, "six")) == (IntegrationMethod::gauss));
            CHECK((StringAccess::GetValue<GridTypeType>(data_set, "seven")) == (GridType::b_spline_grid));
        }

    }  // TEST_SUITE

} // namespace queso::testing
