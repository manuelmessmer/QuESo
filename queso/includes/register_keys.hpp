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

#ifndef REGISTER_KEYS_HPP
#define REGISTER_KEYS_HPP

//// Project includes
#include "queso/includes/keys.hpp"


namespace queso {

// Here we will add all keys.

#if QUESO_BUILD_TESTING
namespace Testing {

    // TestKeys1 :: KeyToSubDict
    QuESo_DEFINE_KEY_SET( TestKeys1, KeyToSubDict, QuESo_KEY_LIST(zero, one, two, three, four) );
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys1, zero);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys1, one);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys1, two);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys1, three);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys1, four);
    QuESo_REGISTER_KEY_SET(TestKeys1, KeyToSubDict,
        QuESo_KEY(TestKeys1::zero),
        QuESo_KEY(TestKeys1::one),
        QuESo_KEY(TestKeys1::two),
        QuESo_KEY(TestKeys1::three),
        QuESo_KEY(TestKeys1::four)
    );

    // TestKeys2 :: KeyToList
    QuESo_DEFINE_KEY_SET( TestKeys2, KeyToList, QuESo_KEY_LIST(zero, one, two, three) );
    QuESo_DEFINE_KEY_TO_LIST(TestKeys2, zero);
    QuESo_DEFINE_KEY_TO_LIST(TestKeys2, one);
    QuESo_DEFINE_KEY_TO_LIST(TestKeys2, two);
    QuESo_DEFINE_KEY_TO_LIST(TestKeys2, three);
    QuESo_REGISTER_KEY_SET(TestKeys2, KeyToList,
        QuESo_KEY(TestKeys2::zero),
        QuESo_KEY(TestKeys2::one),
        QuESo_KEY(TestKeys2::two),
        QuESo_KEY(TestKeys2::three)
    );

    // TestKeys3 :: KeyToValue
    QuESo_DEFINE_KEY_SET( TestKeys3, KeyToValue, QuESo_KEY_LIST(zero, one, two, three, four, five, six, seven) );
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, zero, PointType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, one, Vector3i);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, two, bool);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, three, double);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, four, IndexType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, five, std::string);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, six, IntegrationMethodType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, seven, GridTypeType);
    QuESo_REGISTER_KEY_SET( TestKeys3, KeyToValue,
        QuESo_KEY(TestKeys3::zero),
        QuESo_KEY(TestKeys3::one),
        QuESo_KEY(TestKeys3::two),
        QuESo_KEY(TestKeys3::three),
        QuESo_KEY(TestKeys3::four),
        QuESo_KEY(TestKeys3::five),
        QuESo_KEY(TestKeys3::six),
        QuESo_KEY(TestKeys3::seven)
    );

    // TestKeys4 :: KeyToSubDict
    QuESo_DEFINE_KEY_SET( TestKeys4, KeyToSubDict, QuESo_KEY_LIST(zero, one, two, three) );
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys4, zero);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys4, one);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys4, two);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys4, three);
    QuESo_REGISTER_KEY_SET(TestKeys4, KeyToSubDict,
        QuESo_KEY(TestKeys4::zero),
        QuESo_KEY(TestKeys4::one),
        QuESo_KEY(TestKeys4::two),
        QuESo_KEY(TestKeys4::three)
    );
    // TestKeys4 :: KeyToList
    QuESo_DEFINE_KEY_SET( TestKeys4, KeyToList, QuESo_KEY_LIST(five, six, seven) );
    QuESo_DEFINE_KEY_TO_LIST(TestKeys4, five);
    QuESo_DEFINE_KEY_TO_LIST(TestKeys4, six);
    QuESo_DEFINE_KEY_TO_LIST(TestKeys4, seven);
    QuESo_REGISTER_KEY_SET(TestKeys4, KeyToList,
        QuESo_KEY(TestKeys4::five),
        QuESo_KEY(TestKeys4::six),
        QuESo_KEY(TestKeys4::seven)
    );

    // TestKeys5 :: KeyToSubDict
    QuESo_DEFINE_KEY_SET( TestKeys5, KeyToSubDict, QuESo_KEY_LIST(zero, one) );
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys5, zero);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys5, one);
    QuESo_REGISTER_KEY_SET(TestKeys5, KeyToSubDict,
        QuESo_KEY(TestKeys5::zero),
        QuESo_KEY(TestKeys5::one)
    );
    // TestKeys5 :: KeyToList
    QuESo_DEFINE_KEY_SET( TestKeys5, KeyToList, QuESo_KEY_LIST(five, six) );
    QuESo_DEFINE_KEY_TO_LIST(TestKeys5, five);
    QuESo_DEFINE_KEY_TO_LIST(TestKeys5, six);
    QuESo_REGISTER_KEY_SET(TestKeys5, KeyToList,
        QuESo_KEY(TestKeys5::five),
        QuESo_KEY(TestKeys5::six)
    );
    // TestKeys5 :: KeyToList
    QuESo_DEFINE_KEY_SET( TestKeys5, KeyToValue, QuESo_KEY_LIST(seven, eight, nine, ten) );
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, seven, double);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, eight, IndexType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, nine, PointType );
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, ten, IndexType );
    QuESo_REGISTER_KEY_SET(TestKeys5, KeyToValue,
        QuESo_KEY(TestKeys5::seven),
        QuESo_KEY(TestKeys5::eight),
        QuESo_KEY(TestKeys5::nine),
        QuESo_KEY(TestKeys5::ten)
    );

}
#endif // End QUESO_BUILD_TESTING

} // End namespace queso

#endif // REGISTER_KEYS_HPP