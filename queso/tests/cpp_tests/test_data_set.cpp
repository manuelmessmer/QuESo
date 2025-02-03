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


// This customized tag is needed to get the correct namesapce of GetKeyBaseType().
// template<class T>
// struct FuncForwarder {
//     static auto GetKeyBaseType() -> decltype(queso::Testing::DataSetTestSuite::key::GetKeyBaseType<T>()) {
//     return queso::Testing::DataSetTestSuite::key::GetKeyBaseType<T>();
//     }
// };

// struct FuncForwarder {
//     template<typename TType>
//     static auto GetKeyBaseType() -> decltype(queso::Testing::DataSetTestSuite::key::GetKeyBaseType<TType>()) {
//         return queso::Testing::DataSetTestSuite::key::GetKeyBaseType<TType>();
//     }
// };

BOOST_AUTO_TEST_CASE(TestDataSet) {
    QuESo_INFO << "Testing :: Test DataSet :: ..." << std::endl;

    //typedef decltype(queso::key::GetKeyBaseType<TestKeys19::KeyToDataSet>()) BaseType;

    // DataSet data_set<DataSet::TypeTag<TestKeys19::KeyToDataSet>, FuncForwarder>(DataSet::TypeTag<TestKeys19::KeyToDataSet>{});
    DataSet data_set(DataSet::TypeTag<TestKey::KeyToDataSet>{});
    data_set.SetValue(TestKey::zero, PointType{0.0, 1.0, 2.0});
    // data_set.SetDefault();
    auto a = data_set.GetValue<PointType>(TestKey::zero);
    std::cout << a << std::endl;
    // Data_Set data_set2(Data_Set::TypeTag<Blabla::KeyToDataSet>{});
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso