/*
  ____        ______  _____
 / __ \      |  ____|/ ____|
| |  | |_   _| |__  | (___   ___
| |  | | | | |  __|  \___ \ / _ \
| |__| | |_| | |____ ____) | (_) |
 \___\_\\__,_|______|_____/ \___/
        Quadrature for Embedded Solids

 License:    BSD 4-Clause License
             See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE

 Authors:    Manuel Messmer
*/

#pragma once

//// STL includes
#include <span>

namespace queso {

struct WithoutNormals
{
};

struct WithNormals
{
};

template<typename Mode>
struct TriangleProxy;

template<>
struct TriangleProxy<WithoutNormals>
{
    std::span<const double, 3> P1;
    std::span<const double, 3> P2;
    std::span<const double, 3> P3;
};

template<>
struct TriangleProxy<WithNormals>
{
    std::span<const double, 3> P1;
    std::span<const double, 3> P2;
    std::span<const double, 3> P3;
    std::span<const double, 3> Normal;
};

}// End namespace queso
