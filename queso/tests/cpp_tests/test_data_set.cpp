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

    enum PlainIntEnum {a1, a2, a3};
    enum class ClassIntEnum {a1, a2, a3};

    // Check is_index for
    static_assert( DataSet::is_index_v<IndexType> );
    static_assert( DataSet::is_index_v<unsigned int> );
    static_assert( DataSet::is_index_v<unsigned short> );
    static_assert( DataSet::is_index_v<unsigned long> );
    static_assert( DataSet::is_index_v<int> );
    static_assert( DataSet::is_index_v<short> );
    static_assert( DataSet::is_index_v<long> );

    static_assert( !DataSet::is_index_v<PlainIntEnum> );
    static_assert( !DataSet::is_index_v<ClassIntEnum> );
    static_assert( !DataSet::is_index_v<float> );
    static_assert( !DataSet::is_index_v<double> );
    static_assert( !DataSet::is_index_v<bool> );
    static_assert( !DataSet::is_index_v<char> );
    static_assert( !DataSet::is_index_v<signed char> );
    static_assert( !DataSet::is_index_v<unsigned char> );
    static_assert( !DataSet::is_index_v<std::string> );

    // Check is_unsigned_index
    static_assert( DataSet::is_unsigned_index_v<IndexType> );
    static_assert( DataSet::is_unsigned_index_v<unsigned int> );
    static_assert( DataSet::is_unsigned_index_v<unsigned short> );
    static_assert( DataSet::is_unsigned_index_v<unsigned long> );

    static_assert( !DataSet::is_unsigned_index_v<int> );
    static_assert( !DataSet::is_unsigned_index_v<short> );
    static_assert( !DataSet::is_unsigned_index_v<long> );
    static_assert( !DataSet::is_unsigned_index_v<PlainIntEnum> );
    static_assert( !DataSet::is_unsigned_index_v<ClassIntEnum> );
    static_assert( !DataSet::is_unsigned_index_v<float> );
    static_assert( !DataSet::is_unsigned_index_v<double> );
    static_assert( !DataSet::is_unsigned_index_v<bool> );
    static_assert( !DataSet::is_unsigned_index_v<char> );
    static_assert( !DataSet::is_unsigned_index_v<signed char> );
    static_assert( !DataSet::is_unsigned_index_v<unsigned char> );
    static_assert( !DataSet::is_unsigned_index_v<std::string> );

    // Check is_scoped_int_enum
    static_assert( DataSet::is_scoped_int_enum<ClassIntEnum>::value );

    static_assert( !DataSet::is_scoped_int_enum<PlainIntEnum>::value );
    static_assert( !DataSet::is_scoped_int_enum<int>::value );
    static_assert( !DataSet::is_scoped_int_enum<IndexType>::value );
}

BOOST_AUTO_TEST_CASE(TestDataSet) {
    QuESo_INFO << "Testing :: Test DataSet :: SetValue" << std::endl;

    DataSet data_set(DataSet::TypeTag<key::TestKeys3KeyToValueKeySetInfo>{});
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

    if constexpr (NOTDEBUG) {
        // Wrong keys
        BOOST_REQUIRE_THROW( data_set.SetValue(TestKeys5::seven, 2.0), std::exception );
        BOOST_REQUIRE_THROW( data_set.SetValue(TestKeys5::eight, 2), std::exception );
        BOOST_REQUIRE_THROW( data_set.SetValue(TestKeys5::nine, PointType{1.0, 2.0, 3.0}), std::exception );
        BOOST_REQUIRE_THROW( data_set.SetValue(TestKeys5::ten, 4), std::exception );
    }

    data_set.SetValue("zero", PointType{1.0, 2.0, 3.0});
    data_set.SetValue("zero", Vector3i{1, 2, 3}); // The conversion from Vector3i to PointType is legal.
    BOOST_REQUIRE_THROW( data_set.SetValue("zero", 1.0), std::exception );

    data_set.SetValue("one", Vector3i{1, 2, 3});
    BOOST_REQUIRE_THROW( data_set.SetValue("one", PointType{1.0, 2.0, 3.0}), std::exception );

    data_set.SetValue("two", false);
    BOOST_REQUIRE_THROW( data_set.SetValue("two", 1), std::exception );

    data_set.SetValue("three", 3.0);
    data_set.SetValue("three", 3); // The conversion from int to double is legal.
    data_set.SetValue("three", IndexType(3)); // The conversion from int to double is legal.
    BOOST_REQUIRE_THROW( data_set.SetValue("three", true), std::exception );

    data_set.SetValue("four", (int)1);
    data_set.SetValue("four", (unsigned int)1 );
    data_set.SetValue("four", IndexType(1));
    data_set.SetValue("four", (long)1 );
    data_set.SetValue("four", (short)1 );
    BOOST_REQUIRE_THROW( data_set.SetValue("four", -1), std::exception );
    BOOST_REQUIRE_THROW( data_set.SetValue("four", 1.0), std::exception );

    data_set.SetValue("five", std::string("hallo"));
    BOOST_REQUIRE_THROW( data_set.SetValue("five", 1), std::exception );

    data_set.SetValue("six", IntegrationMethod::gauss);
    BOOST_REQUIRE_THROW( data_set.SetValue("six", 1), std::exception );

    data_set.SetValue("seven", GridType::b_spline_grid);
    BOOST_REQUIRE_THROW( data_set.SetValue("seven", 1), std::exception );

    BOOST_REQUIRE_THROW( data_set.SetValue("wrong_key", 1), std::exception );
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso