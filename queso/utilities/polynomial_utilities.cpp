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

#include "queso/utilities/polynomial_utilities.h"
#include "queso/includes/define.hpp"

namespace queso {

template<typename T, typename... Args>
auto HelperMakeVector(Args&&... args) {
  std::vector<T> vec;
  vec.reserve(sizeof...(Args));
  (vec.emplace_back(std::forward<Args>(args)), ...);
  return vec;
}

const std::vector<Unique<Polynomial::FuncBase>> Polynomial::mLegendePolynomials =
    HelperMakeVector<Unique<Polynomial::FuncBase>>(
        MakeUnique<FxP0>(),
        MakeUnique<FxP1>(),
        MakeUnique<FxP2>(),
        MakeUnique<FxP3>(),
        MakeUnique<FxP4>(),
        MakeUnique<FxP5>(),
        MakeUnique<FxP6>(),
        MakeUnique<FxP7>(),
        MakeUnique<FxP8>() );

} // End namespace queso

