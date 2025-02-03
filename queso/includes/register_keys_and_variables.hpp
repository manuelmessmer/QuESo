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
#include "queso/includes/variables.hpp"

namespace queso {

// Here we will register all Keys and Variables used in QuESo
QuESo_REGISTER_KEY_SET_1(DummyKey, DataSet, QuESo_LIST(zero, one) );

QuESo_REGISTER_KEY_SET_1(TestKey, DataSet, QuESo_LIST(zero, one) );
QuESo_REGISTER_DATASET_VARIABLES(TestKey, 
    QuESo_VARIABLE(TestKey::zero, PointType),
    QuESo_VARIABLE(TestKey::one, int),
);

} // End namespace queso

#endif // REGISTER_KEYS_AND_VARIABLES