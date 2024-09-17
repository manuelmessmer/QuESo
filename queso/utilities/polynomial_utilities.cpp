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

// Project includes
#include "queso/utilities/polynomial_utilities.h"

namespace queso {


double Polynomial::f_x(double x, int order, double a, double b){
        return mLegendePolynomials[order]->f_x(x, a, b);
    }

double Polynomial::f_x_int(double x, int order, double a, double b){
    return mLegendePolynomials[order]->f_x_int(x, a, b);
}

template<typename T, typename... Args>
auto HelperMakeVector(Args&&... args) {
  std::vector<T> vec;
  vec.reserve(sizeof...(Args));
  (vec.emplace_back(std::forward<Args>(args)), ...);
  return vec;
}

const std::vector<Unique<Polynomial::FuncBase>> Polynomial::mLegendePolynomials =
    HelperMakeVector<Unique<Polynomial::FuncBase>>(
        MakeUnique<Polynomial::FxP0>(),
        MakeUnique<Polynomial::FxP1>(),
        MakeUnique<Polynomial::FxP2>(),
        MakeUnique<Polynomial::FxP3>(),
        MakeUnique<Polynomial::FxP4>(),
        MakeUnique<Polynomial::FxP5>(),
        MakeUnique<Polynomial::FxP6>(),
        MakeUnique<Polynomial::FxP7>(),
        MakeUnique<Polynomial::FxP8>() );

} // End namespace queso

