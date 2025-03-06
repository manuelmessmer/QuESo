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

// Create value type lists
QuESo_CREATE_VALUE_TYPE_LIST(MainValueTypeList, PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType, GridTypeType);

/// Create type tags.
// We could create different KeySetToValues's here. However, for now, one ("Main") seems sufficient.
QuESo_CREATE_KEY_SET_TO_VALUE_TYPE_TAG(MainValuesTypeTag, MainValueTypeList);
QuESo_CREATE_KEY_SET_TO_OBJECT_TYPE_TAG(ListTypeTag);
QuESo_CREATE_KEY_SET_TO_OBJECT_TYPE_TAG(SubDictTypeTag);

// Here, we will add all global keys.

#if QUESO_BUILD_TESTING
namespace Testing {

    // TestKeys1 :: SubDictTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys1, SubDictTypeTag, QuESo_KEY_LIST(zero, one, two, three, four) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys1, zero, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys1, one, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys1, two, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys1, three, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys1, four, SubDictTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys1, SubDictTypeTag,
        QuESo_KEY(TestKeys1::zero),
        QuESo_KEY(TestKeys1::one),
        QuESo_KEY(TestKeys1::two),
        QuESo_KEY(TestKeys1::three),
        QuESo_KEY(TestKeys1::four)
    );

    // TestKeys2 :: ListTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys2, ListTypeTag, QuESo_KEY_LIST(zero, one, two, three) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys2, zero, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys2, one, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys2, two, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys2, three, ListTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys2, ListTypeTag,
        QuESo_KEY(TestKeys2::zero),
        QuESo_KEY(TestKeys2::one),
        QuESo_KEY(TestKeys2::two),
        QuESo_KEY(TestKeys2::three)
    );

    // TestKeys3 :: MainValuesTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys3, MainValuesTypeTag, QuESo_KEY_LIST(zero, one, two, three, four, five, six, seven) );
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, zero, MainValuesTypeTag, PointType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, one, MainValuesTypeTag, Vector3i);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, two, MainValuesTypeTag, bool);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, three, MainValuesTypeTag, double);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, four, MainValuesTypeTag, IndexType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, five, MainValuesTypeTag, std::string);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, six, MainValuesTypeTag, IntegrationMethodType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, seven, MainValuesTypeTag, GridTypeType);
    QuESo_REGISTER_KEY_SET( TestKeys3, MainValuesTypeTag,
        QuESo_KEY(TestKeys3::zero),
        QuESo_KEY(TestKeys3::one),
        QuESo_KEY(TestKeys3::two),
        QuESo_KEY(TestKeys3::three),
        QuESo_KEY(TestKeys3::four),
        QuESo_KEY(TestKeys3::five),
        QuESo_KEY(TestKeys3::six),
        QuESo_KEY(TestKeys3::seven)
    );

    // TestKeys4 :: SubDictTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys4, SubDictTypeTag, QuESo_KEY_LIST(zero, one, two, three) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, zero, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, one, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, two, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, three, SubDictTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys4, SubDictTypeTag,
        QuESo_KEY(TestKeys4::zero),
        QuESo_KEY(TestKeys4::one),
        QuESo_KEY(TestKeys4::two),
        QuESo_KEY(TestKeys4::three)
    );
    // TestKeys4 :: ListTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys4, ListTypeTag, QuESo_KEY_LIST(five, six, seven) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, five, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, six, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, seven, ListTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys4, ListTypeTag,
        QuESo_KEY(TestKeys4::five),
        QuESo_KEY(TestKeys4::six),
        QuESo_KEY(TestKeys4::seven)
    );

    // TestKeys5 :: SubDictTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys5, SubDictTypeTag, QuESo_KEY_LIST(zero, one) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys5, zero, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys5, one, SubDictTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys5, SubDictTypeTag,
        QuESo_KEY(TestKeys5::zero),
        QuESo_KEY(TestKeys5::one)
    );
    // TestKeys5 :: ListTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys5, ListTypeTag, QuESo_KEY_LIST(five, six) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys5, five, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys5, six, ListTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys5, ListTypeTag,
        QuESo_KEY(TestKeys5::five),
        QuESo_KEY(TestKeys5::six)
    );
    // TestKeys5 :: ListTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys5, MainValuesTypeTag, QuESo_KEY_LIST(seven, eight, nine, ten) );
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, seven, MainValuesTypeTag, double);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, eight, MainValuesTypeTag, IndexType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, nine, MainValuesTypeTag, PointType );
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, ten, MainValuesTypeTag, IndexType );
    QuESo_REGISTER_KEY_SET(TestKeys5, MainValuesTypeTag,
        QuESo_KEY(TestKeys5::seven),
        QuESo_KEY(TestKeys5::eight),
        QuESo_KEY(TestKeys5::nine),
        QuESo_KEY(TestKeys5::ten)
    );

}
#endif // End QUESO_BUILD_TESTING

} // End namespace queso

#endif // REGISTER_KEYS_HPP