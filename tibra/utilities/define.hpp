// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef DEFINE_INCLUDE_HPP
#define DEFINE_INCLUDE_HPP

//// STL includes
#include <limits>
#include <memory>

namespace tibra {
///@name TIBRA GLOBAL VARIABLES
///@{

// Large negative number
constexpr double LOWESTD = std::numeric_limits<double>::lowest();
// Basically zero, very close to zero, but positive.
constexpr double MIND = std::numeric_limits<double>::min();
// Large positive number
constexpr double MAXD = std::numeric_limits<double>::max();
// Infinity (even larger than max)
constexpr double INFINITYD = std::numeric_limits<double>::infinity();
// Is NaN
constexpr double QUIETNAND = std::numeric_limits<double>::quiet_NaN();

// Tolerances
constexpr double EPS0 = 1e-7;
constexpr double EPS1 = 1e-8;
constexpr double EPS2 = 1e-10;
constexpr double EPS3 = 1e-12;
constexpr double EPS4 = 1e-14;

///@}
///@name TIBRA POINTER DEFINITIONS
///@{

// Pointer Definitions
template <typename T>
using Shared = std::shared_ptr<T>;

template <typename T, typename... Args>
auto MakeShared(Args&&... args) -> decltype(std::make_shared<T>(std::forward<Args>(args)...)) {
  return std::make_shared<T>(std::forward<Args>(args)...);
}

template <typename T>
using Unique = std::unique_ptr<T>;

template <typename T, typename... Args>
auto MakeUnique(Args&&... args) -> decltype(std::make_unique<T>(std::forward<Args>(args)...)) {
  return std::make_unique<T>(std::forward<Args>(args)...);
}

///@}
///@name TIBRA GLOBAL OPERATIONS
///@{

namespace op {

}

///@}
} // End namespace tibra

#endif // DEFINE_INCLUDE_HPP