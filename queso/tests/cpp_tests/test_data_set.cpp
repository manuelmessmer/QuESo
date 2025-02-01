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
#include "queso/containers/data_set.hpp"

namespace queso {

namespace Testing {

BOOST_AUTO_TEST_SUITE( DataSetTestSuite )

QuESo_REGISTER_KEY_SET_1(TestKeys19, DataSet, QuESo_LIST(zero, one) );

// This customized tag is needed to get the correct namesapce of GetKeyBaseType().
template<class T>
struct TypeTagTester {
    static auto GetKeyBaseType() -> decltype(queso::Testing::DataSetTestSuite::key::GetKeyBaseType<T>()) {
    return queso::Testing::DataSetTestSuite::key::GetKeyBaseType<T>();
    }
};

BOOST_AUTO_TEST_CASE(TestDataSet) {
    QuESo_INFO << "Testing :: Test DataSet :: ..." << std::endl;

    //typedef decltype(queso::key::GetKeyBaseType<TestKeys19::KeyToDataSet>()) BaseType;

    Data_Set data_set(TypeTagTester<TestKeys19::KeyToDataSet>{});
    // data_set.SetDefault();

    // Data_Set data_set2(Data_Set::TypeTag<Blabla::KeyToDataSet>{});
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso