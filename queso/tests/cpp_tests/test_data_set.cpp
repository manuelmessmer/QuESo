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
#include "queso/containers/data_set.hpp"

namespace queso {

namespace Testing {

BOOST_AUTO_TEST_SUITE( DataSetTestSuite )

BOOST_AUTO_TEST_CASE(TestDataSetTypeTraits) {
    QuESo_INFO << "Testing :: Test DataSet :: Type Traits" << std::endl;

    using DataSetType = DataSet<queso::key::MainValuesTypeTag>;

    enum PlainIntEnum {a1, a2, a3};
    enum class ClassIntEnum {a1, a2, a3};

    // Check is_index for
    static_assert( DataSetType::is_index_v<IndexType> );
    static_assert( DataSetType::is_index_v<unsigned int> );
    static_assert( DataSetType::is_index_v<unsigned short> );
    static_assert( DataSetType::is_index_v<unsigned long> );
    static_assert( DataSetType::is_index_v<int> );
    static_assert( DataSetType::is_index_v<short> );
    static_assert( DataSetType::is_index_v<long> );

    static_assert( !DataSetType::is_index_v<PlainIntEnum> );
    static_assert( !DataSetType::is_index_v<ClassIntEnum> );
    static_assert( !DataSetType::is_index_v<float> );
    static_assert( !DataSetType::is_index_v<double> );
    static_assert( !DataSetType::is_index_v<bool> );
    static_assert( !DataSetType::is_index_v<char> );
    static_assert( !DataSetType::is_index_v<signed char> );
    static_assert( !DataSetType::is_index_v<unsigned char> );
    static_assert( !DataSetType::is_index_v<std::string> );

    // Check is_unsigned_index
    static_assert( DataSetType::is_unsigned_index_v<IndexType> );
    static_assert( DataSetType::is_unsigned_index_v<unsigned int> );
    static_assert( DataSetType::is_unsigned_index_v<unsigned short> );
    static_assert( DataSetType::is_unsigned_index_v<unsigned long> );

    static_assert( !DataSetType::is_unsigned_index_v<int> );
    static_assert( !DataSetType::is_unsigned_index_v<short> );
    static_assert( !DataSetType::is_unsigned_index_v<long> );
    static_assert( !DataSetType::is_unsigned_index_v<PlainIntEnum> );
    static_assert( !DataSetType::is_unsigned_index_v<ClassIntEnum> );
    static_assert( !DataSetType::is_unsigned_index_v<float> );
    static_assert( !DataSetType::is_unsigned_index_v<double> );
    static_assert( !DataSetType::is_unsigned_index_v<bool> );
    static_assert( !DataSetType::is_unsigned_index_v<char> );
    static_assert( !DataSetType::is_unsigned_index_v<signed char> );
    static_assert( !DataSetType::is_unsigned_index_v<unsigned char> );
    static_assert( !DataSetType::is_unsigned_index_v<std::string> );

    // Check is_scoped_int_enum
    static_assert( DataSetType::is_key_v<decltype(TestKeys3::zero)> );

    static_assert( !DataSetType::is_key_v<ClassIntEnum> );
    static_assert( !DataSetType::is_key_v<PlainIntEnum> );
    static_assert( !DataSetType::is_key_v<int> );
    static_assert( !DataSetType::is_key_v<IndexType> );
}

BOOST_AUTO_TEST_CASE(TestDataSetValue) {
    QuESo_INFO << "Testing :: Test DataSet :: SetValue" << std::endl;

    using DataSetType = DataSet<queso::key::MainValuesTypeTag>;
    DataSetType data_set(DataSetType::KeySetInfoTypeTag<key::detail::TestKeys3MainValuesTypeTagKeySetInfo>{});

    data_set.SetValue(TestKeys3::zero, PointType{1.0, 2.0, 3.0});
    data_set.SetValue(TestKeys3::one, Vector3i{1, 2, 3});
    data_set.SetValue(TestKeys3::two, true);
    data_set.SetValue(TestKeys3::three, 2.0);
    data_set.SetValue(TestKeys3::four, (unsigned int)1 );
    data_set.SetValue(TestKeys3::four, (unsigned short)1 );
    data_set.SetValue(TestKeys3::four, (unsigned long)1 );
    data_set.SetValue(TestKeys3::four, IndexType(1));
    data_set.SetValue(TestKeys3::five, std::string("Hallo"));
    data_set.SetValue(TestKeys3::six, IntegrationMethod::gauss);
    data_set.SetValue(TestKeys3::seven, GridType::b_spline_grid);

    if constexpr (!NOTDEBUG) {
        // Wrong keys
        BOOST_REQUIRE_THROW( data_set.SetValue(TestKeys5::seven, 2.0), std::exception );
        BOOST_REQUIRE_THROW( data_set.SetValue(TestKeys5::eight, IndexType(2)), std::exception );
        BOOST_REQUIRE_THROW( data_set.SetValue(TestKeys5::nine, PointType{1.0, 2.0, 3.0}), std::exception );
        BOOST_REQUIRE_THROW( data_set.SetValue(TestKeys5::ten, 4ul), std::exception );
    }

    using StringAccess = DataSetStringAccess<DataSetType>;

    StringAccess::SetValue(data_set, "zero", PointType{1.0, 2.0, 3.0});
    StringAccess::SetValue(data_set, "zero", Vector3i{1, 2, 3}); // The conversion from Vector3i to PointType is legal.
    BOOST_REQUIRE_THROW( StringAccess::SetValue(data_set, "zero", 1.0), std::exception );

    StringAccess::SetValue( data_set, "one", Vector3i{1, 2, 3});
    BOOST_REQUIRE_THROW( StringAccess::SetValue(data_set, "one", PointType{1.0, 2.0, 3.0}), std::exception );

    StringAccess::SetValue(data_set, "two", false);
    BOOST_REQUIRE_THROW( StringAccess::SetValue(data_set, "two", 1), std::exception );

    StringAccess::SetValue(data_set, "three", 3.0);
    StringAccess::SetValue(data_set, "three", 3); // The conversion from int to double is legal.
    StringAccess::SetValue(data_set, "three", IndexType(3)); // The conversion from int to double is legal.
    BOOST_REQUIRE_THROW( StringAccess::SetValue(data_set, "three", true), std::exception );

    StringAccess::SetValue(data_set, "four", (int)1);
    StringAccess::SetValue(data_set, "four", (unsigned int)1 );
    StringAccess::SetValue(data_set, "four", IndexType(1));
    StringAccess::SetValue(data_set, "four", (long)1 );
    StringAccess::SetValue(data_set, "four", (short)1 );
    BOOST_REQUIRE_THROW( StringAccess::SetValue(data_set, "four", -1), std::exception );
    BOOST_REQUIRE_THROW( StringAccess::SetValue(data_set, "four", 1.0), std::exception );

    StringAccess::SetValue(data_set, "five", std::string("hallo"));
    BOOST_REQUIRE_THROW( StringAccess::SetValue(data_set, "five", 1), std::exception );

    StringAccess::SetValue(data_set, "six", IntegrationMethod::gauss);
    BOOST_REQUIRE_THROW( StringAccess::SetValue(data_set, "six", Vector3i{0, 0, 0}), std::exception );

    StringAccess::SetValue(data_set, "seven", GridType::b_spline_grid);
    BOOST_REQUIRE_THROW( StringAccess::SetValue(data_set, "seven", 1), std::exception );

    BOOST_REQUIRE_THROW( StringAccess::SetValue(data_set, "wrong_key", 1), std::exception );
}

BOOST_AUTO_TEST_CASE(TestDataSetGetValueKey) {
    QuESo_INFO << "Testing :: Test DataSet :: Set/GetValue (Key)" << std::endl;

    using DataSetType = DataSet<queso::key::MainValuesTypeTag>;
    DataSetType data_set(DataSetType::KeySetInfoTypeTag<key::detail::TestKeys3MainValuesTypeTagKeySetInfo>{});

    // Values are not set.
    QuESo_CHECK( !data_set.IsSet(TestKeys3::zero) );
    QuESo_CHECK( !data_set.IsSet(TestKeys3::one) );
    QuESo_CHECK( !data_set.IsSet(TestKeys3::two) );
    QuESo_CHECK( !data_set.IsSet(TestKeys3::three) );
    QuESo_CHECK( !data_set.IsSet(TestKeys3::four) );
    QuESo_CHECK( !data_set.IsSet(TestKeys3::five) );
    QuESo_CHECK( !data_set.IsSet(TestKeys3::six) );
    QuESo_CHECK( !data_set.IsSet(TestKeys3::seven) );

    // Values are not set.
    BOOST_REQUIRE_THROW( data_set.GetValue<PointType>(TestKeys3::zero), std::exception);
    BOOST_REQUIRE_THROW( data_set.GetValue<Vector3i>(TestKeys3::one), std::exception);
    BOOST_REQUIRE_THROW( data_set.GetValue<bool>(TestKeys3::two), std::exception);
    BOOST_REQUIRE_THROW( data_set.GetValue<double>(TestKeys3::three), std::exception);
    BOOST_REQUIRE_THROW( data_set.GetValue<IndexType>(TestKeys3::four), std::exception);
    BOOST_REQUIRE_THROW( data_set.GetValue<std::string>(TestKeys3::five), std::exception);
    BOOST_REQUIRE_THROW( data_set.GetValue<IntegrationMethodType>(TestKeys3::six), std::exception);
    BOOST_REQUIRE_THROW( data_set.GetValue<GridTypeType>(TestKeys3::seven), std::exception);

    // Values are not set.
    if constexpr(!NOTDEBUG) {
        BOOST_REQUIRE_THROW( data_set.GetValueFast<PointType>(TestKeys3::zero), std::exception);
        BOOST_REQUIRE_THROW( data_set.GetValueFast<Vector3i>(TestKeys3::one), std::exception);
        BOOST_REQUIRE_THROW( data_set.GetValueFast<bool>(TestKeys3::two), std::exception);
        BOOST_REQUIRE_THROW( data_set.GetValueFast<double>(TestKeys3::three), std::exception);
        BOOST_REQUIRE_THROW( data_set.GetValueFast<IndexType>(TestKeys3::four), std::exception);
        BOOST_REQUIRE_THROW( data_set.GetValueFast<std::string>(TestKeys3::five), std::exception);
        BOOST_REQUIRE_THROW( data_set.GetValueFast<IntegrationMethodType>(TestKeys3::six), std::exception);
        BOOST_REQUIRE_THROW( data_set.GetValueFast<GridTypeType>(TestKeys3::seven), std::exception);
    }

    // Set values
    data_set.SetValue(TestKeys3::zero, PointType{1.0, 2.0, 3.0});
    data_set.SetValue(TestKeys3::one, Vector3i{1, 2, 3});
    data_set.SetValue(TestKeys3::two, true);
    data_set.SetValue(TestKeys3::three, 2.0);
    data_set.SetValue(TestKeys3::four, IndexType(5));
    data_set.SetValue(TestKeys3::five, std::string("Hallo"));
    data_set.SetValue(TestKeys3::six, IntegrationMethod::gauss);
    data_set.SetValue(TestKeys3::seven, GridType::b_spline_grid);

    // Values are set.
    QuESo_CHECK( data_set.IsSet(TestKeys3::zero) );
    QuESo_CHECK( data_set.IsSet(TestKeys3::one) );
    QuESo_CHECK( data_set.IsSet(TestKeys3::two) );
    QuESo_CHECK( data_set.IsSet(TestKeys3::three) );
    QuESo_CHECK( data_set.IsSet(TestKeys3::four) );
    QuESo_CHECK( data_set.IsSet(TestKeys3::five) );
    QuESo_CHECK( data_set.IsSet(TestKeys3::six) );
    QuESo_CHECK( data_set.IsSet(TestKeys3::seven) );

    // Wrong keys
    if constexpr (!NOTDEBUG) {
        BOOST_REQUIRE_THROW( data_set.IsSet(TestKeys5::seven), std::exception );

        BOOST_REQUIRE_THROW( data_set.GetValue<double>(TestKeys5::seven), std::exception );
        BOOST_REQUIRE_THROW( data_set.GetValue<IndexType>(TestKeys5::eight), std::exception );
        BOOST_REQUIRE_THROW( data_set.GetValue<PointType>(TestKeys5::nine), std::exception );
        BOOST_REQUIRE_THROW( data_set.GetValue<IndexType>(TestKeys5::ten), std::exception );

        BOOST_REQUIRE_THROW( data_set.GetValueFast<double>(TestKeys5::seven), std::exception );
        BOOST_REQUIRE_THROW( data_set.GetValueFast<IndexType>(TestKeys5::eight), std::exception );
        BOOST_REQUIRE_THROW( data_set.GetValueFast<PointType>(TestKeys5::nine), std::exception );
        BOOST_REQUIRE_THROW( data_set.GetValueFast<IndexType>(TestKeys5::ten), std::exception );
    }

    // Get values
    QuESo_CHECK_POINT_NEAR( data_set.GetValue<PointType>(TestKeys3::zero), PointType({1.0, 2.0, 3.0}), 1e-12);
    QuESo_CHECK_Vector3i_EQUAL( data_set.GetValue<Vector3i>(TestKeys3::one), Vector3i({1, 2, 3}) );
    QuESo_CHECK_EQUAL( data_set.GetValue<bool>(TestKeys3::two), true );
    QuESo_CHECK_NEAR( data_set.GetValue<double>(TestKeys3::three), 2.0, 1e-12);
    QuESo_CHECK_EQUAL( data_set.GetValue<IndexType>(TestKeys3::four), 5 );
    QuESo_CHECK_EQUAL( data_set.GetValue<std::string>(TestKeys3::five), std::string("Hallo") );
    QuESo_CHECK_EQUAL( data_set.GetValue<IntegrationMethodType>(TestKeys3::six), IntegrationMethod::gauss );
    QuESo_CHECK_EQUAL( data_set.GetValue<GridTypeType>(TestKeys3::seven), GridType::b_spline_grid );
}


BOOST_AUTO_TEST_CASE(TestDataSetGetValueKeyName) {
    QuESo_INFO << "Testing :: Test DataSet :: Set/GetValue (KeyName)" << std::endl;

    using DataSetType = DataSet<queso::key::MainValuesTypeTag>;
    DataSetType data_set(DataSetType::KeySetInfoTypeTag<key::detail::TestKeys3MainValuesTypeTagKeySetInfo>{});

    using StringAccess = DataSetStringAccess<DataSetType>;

    // Has keys.
    QuESo_CHECK( StringAccess::Has(data_set, "zero") );
    QuESo_CHECK( StringAccess::Has(data_set, "one") );
    QuESo_CHECK( StringAccess::Has(data_set, "two") );
    QuESo_CHECK( StringAccess::Has(data_set, "three") );
    QuESo_CHECK( StringAccess::Has(data_set, "four") );
    QuESo_CHECK( StringAccess::Has(data_set, "five") );
    QuESo_CHECK( StringAccess::Has(data_set, "six") );
    QuESo_CHECK( StringAccess::Has(data_set, "seven") );

    // Does not have key.
    QuESo_CHECK( !StringAccess::Has(data_set, "wrong_key") );

    // Values are not set.
    QuESo_CHECK( !StringAccess::IsSet(data_set, "zero") );
    QuESo_CHECK( !StringAccess::IsSet(data_set, "one") );
    QuESo_CHECK( !StringAccess::IsSet(data_set, "two") );
    QuESo_CHECK( !StringAccess::IsSet(data_set, "three") );
    QuESo_CHECK( !StringAccess::IsSet(data_set, "four") );
    QuESo_CHECK( !StringAccess::IsSet(data_set, "five") );
    QuESo_CHECK( !StringAccess::IsSet(data_set, "six") );
    QuESo_CHECK( !StringAccess::IsSet(data_set, "seven") );

    // Values are not set.
    BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(data_set, "zero"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<Vector3i>(data_set, "one"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<bool>(data_set, "two"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(data_set, "three"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(data_set, "four"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<std::string>(data_set, "five"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<IntegrationMethodType>(data_set, "six"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<GridTypeType>(data_set, "seven"), std::exception);

    // Set values
    StringAccess::SetValue(data_set, "zero", PointType{1.0, 2.0, 3.0});
    StringAccess::SetValue(data_set, "one", Vector3i{1, 2, 3});
    StringAccess::SetValue(data_set, "two", true);
    StringAccess::SetValue(data_set, "three", 2.0);
    StringAccess::SetValue(data_set, "four", IndexType(5));
    StringAccess::SetValue(data_set, "five", std::string("Hallo"));
    StringAccess::SetValue(data_set, "six", IntegrationMethod::gauss);
    StringAccess::SetValue(data_set, "seven", GridType::b_spline_grid);

    // Values are set.
    QuESo_CHECK( StringAccess::IsSet(data_set, "zero") );
    QuESo_CHECK( StringAccess::IsSet(data_set, "one") );
    QuESo_CHECK( StringAccess::IsSet(data_set, "two") );
    QuESo_CHECK( StringAccess::IsSet(data_set, "three") );
    QuESo_CHECK( StringAccess::IsSet(data_set, "four") );
    QuESo_CHECK( StringAccess::IsSet(data_set, "five") );
    QuESo_CHECK( StringAccess::IsSet(data_set, "six") );
    QuESo_CHECK( StringAccess::IsSet(data_set, "seven") );

    // Wrong key
    BOOST_REQUIRE_THROW( StringAccess::IsSet(data_set, "zero_"), std::exception);

    // Wrong keys
    BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(data_set, "zero_"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<Vector3i>(data_set, "one_"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<bool>(data_set, "two_"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(data_set, "three_"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(data_set, "four_"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<std::string>(data_set, "five_"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<IntegrationMethodType>(data_set, "six_"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<GridTypeType>(data_set, "seven_"), std::exception);

    // Wrong variable types
    BOOST_REQUIRE_THROW( StringAccess::GetValue<Vector3i>(data_set, "zero"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<bool>(data_set, "one"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(data_set, "two"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(data_set, "three"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(data_set, "four"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<GridTypeType>(data_set, "five"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<std::string>(data_set, "six"), std::exception);
    BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(data_set, "seven"), std::exception);

    // Get values
    QuESo_CHECK_POINT_NEAR( StringAccess::GetValue<PointType>(data_set, "zero"), PointType({1.0, 2.0, 3.0}), 1e-12);
    QuESo_CHECK_Vector3i_EQUAL( StringAccess::GetValue<Vector3i>(data_set, "one"), Vector3i({1, 2, 3}) );
    QuESo_CHECK_EQUAL( StringAccess::GetValue<bool>(data_set, "two"), true );
    QuESo_CHECK_NEAR( StringAccess::GetValue<double>(data_set, "three"), 2.0, 1e-12);
    QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(data_set, "four"), 5 );
    QuESo_CHECK_EQUAL( StringAccess::GetValue<std::string>(data_set, "five"), std::string("Hallo") );
    QuESo_CHECK_EQUAL( StringAccess::GetValue<IntegrationMethodType>(data_set, "six"), IntegrationMethod::gauss );
    QuESo_CHECK_EQUAL( StringAccess::GetValue<GridTypeType>(data_set, "seven"), GridType::b_spline_grid );
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso