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

#ifndef REGISTER_KEYS_AND_VARIABLES
#define REGISTER_KEYS_AND_VARIABLES

//// Project includes
#include "queso/includes/keys.hpp"


namespace queso {

QuESo_DEFINE_KEY_SET(DummyKey, KeyToValue, QuESo_LIST(zero, one) );
QuESo_DEFINE_KEY_SET(DummyKey, KeyToSubDict, QuESo_LIST(two, three) );
QuESo_DEFINE_KEY_SET(DummyKey, KeyToList, QuESo_LIST(four, five) );

QuESo_DEFINE_KEY_TO_VALUE(DummyKey, zero, PointType);
QuESo_DEFINE_KEY_TO_VALUE(DummyKey, one, int);

QuESo_DEFINE_KEY_TO_SUBDICT(DummyKey, two);
QuESo_DEFINE_KEY_TO_SUBDICT(DummyKey, three);

QuESo_DEFINE_KEY_TO_LIST(DummyKey, four);
QuESo_DEFINE_KEY_TO_LIST(DummyKey, five);

namespace Testing {
    
    QuESo_DEFINE_KEY_SET( TestKeys1, KeyToSubDict, QuESo_LIST(zero, one, two, three, four) );
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys1, zero);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys1, one);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys1, two);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys1, three);
    QuESo_DEFINE_KEY_TO_SUBDICT(TestKeys1, four);

    QuESo_DEFINE_KEY_SET( TestKeys2, KeyToList, QuESo_LIST(zero, one, two, three) );
    QuESo_DEFINE_KEY_TO_LIST(TestKeys2, zero);
    QuESo_DEFINE_KEY_TO_LIST(TestKeys2, one);
    QuESo_DEFINE_KEY_TO_LIST(TestKeys2, two);
    QuESo_DEFINE_KEY_TO_LIST(TestKeys2, three);

    QuESo_DEFINE_KEY_SET( TestKeys3, KeyToValue, QuESo_LIST(zero, one) ); 
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, zero, PointType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, one, int);

}

} // End namespace queso

#endif // REGISTER_KEYS_AND_VARIABLES