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

//// Project includes
#include "queso/includes/register_keys.hpp"


namespace queso {

// cpp
QuESo_REGISTER_KEY_SET(DummyKey, KeyToValue, 
    QuESo_KEY(DummyKey::zero),
    QuESo_KEY(DummyKey::one)
)

QuESo_REGISTER_KEY_SET(DummyKey, KeyToSubDict, 
    QuESo_KEY(DummyKey::two),
    QuESo_KEY(DummyKey::three)
)

QuESo_REGISTER_KEY_SET(DummyKey, KeyToList, 
    QuESo_KEY(DummyKey::four),
    QuESo_KEY(DummyKey::five)
)


namespace Testing {
    
    QuESo_CREATE_KEY(TestKeys1, zero);
    QuESo_CREATE_KEY(TestKeys1, one);
    QuESo_CREATE_KEY(TestKeys1, two);
    QuESo_CREATE_KEY(TestKeys1, three);
    QuESo_CREATE_KEY(TestKeys1, four);
    QuESo_REGISTER_KEY_SET( TestKeys1, KeyToSubDict, 
        QuESo_KEY(TestKeys1::zero),
        QuESo_KEY(TestKeys1::one),
        QuESo_KEY(TestKeys1::two),
        QuESo_KEY(TestKeys1::three),
        QuESo_KEY(TestKeys1::four)
    )

    QuESo_CREATE_KEY(TestKeys2, zero);
    QuESo_CREATE_KEY(TestKeys2, one);
    QuESo_CREATE_KEY(TestKeys2, two);
    QuESo_CREATE_KEY(TestKeys2, three);
    QuESo_REGISTER_KEY_SET( TestKeys2, KeyToList, 
        QuESo_KEY(TestKeys2::zero),
        QuESo_KEY(TestKeys2::one),
        QuESo_KEY(TestKeys2::two),
        QuESo_KEY(TestKeys2::three)
    )

    QuESo_CREATE_KEY(TestKeys3, zero);
    QuESo_CREATE_KEY(TestKeys3, one);
    QuESo_REGISTER_KEY_SET( TestKeys3, KeyToValue, 
        QuESo_KEY(TestKeys3::zero),
        QuESo_KEY(TestKeys3::one)
    )

}

} // End namespace queso
