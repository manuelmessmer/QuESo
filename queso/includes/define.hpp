// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef DEFINE_INCLUDE_HPP
#define DEFINE_INCLUDE_HPP

//// STL includes
#include <limits>
#include <memory>
//// Project includes
#include "includes/logger.hpp"
#include "includes/exception.hpp"
#include "includes/timer.hpp"
#include "containers/vector3.hpp"

namespace queso {

///@name QuESo GLOBAL TYPE DEFINITIONS
///@{

typedef std::size_t  SizeType;
typedef std::size_t  IndexType;

typedef Vector3<double> PointType;
typedef Vector3<double> Vector3d;
typedef Vector3<IndexType> Vector3i;
typedef std::pair<PointType, PointType> BoundingBoxType;

///@}
///@name QuESo GLOBAL VARIABLES
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
constexpr double SNAPTOL = 1e-12;
constexpr double ZEROTOL = 1e-14;

inline double RelativeSnapTolerance(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance = SNAPTOL){
    const auto delta = (rUpperBound - rLowerBound);
    return std::max( std::max(delta[0], std::max(delta[1], delta[2]))*Tolerance, Tolerance);
}

inline double RelativeSnapTolerance(const PointType& rDelta, double Tolerance = SNAPTOL){
    return std::max( std::max(rDelta[0], std::max(rDelta[1], rDelta[2]))*Tolerance, Tolerance);
}

// Tolerances for Testing
constexpr double EPS0 = 1e-7;
constexpr double EPS1 = 1e-8;
constexpr double EPS2 = 1e-10;
constexpr double EPS3 = 1e-12;
constexpr double EPS4 = 1e-14;

// Enum's
enum IntegrationMethod {Gauss, Gauss_Reduced1, Gauss_Reduced2, GGQ_Optimal, GGQ_Reduced1, GGQ_Reduced2};
typedef enum IntegrationMethod IntegrationMethodType;

enum IntersectionStatus {Inside, Outside, Trimmed};
typedef IntersectionStatus IntersectionStatusType;

// QuESo Factories
inline BoundingBoxType MakeBox( PointType rL, PointType rR ){
    return std::make_pair(rL, rR);
}

///@}
///@name QuESo POINTER DEFINITIONS
///@{

// Shared Ptr
template <typename T>
using Shared = std::shared_ptr<T>;

template <typename T, typename... Args>
auto MakeShared(Args&&... args) -> decltype(std::make_shared<T>(std::forward<Args>(args)...)) {
    return std::make_shared<T>(std::forward<Args>(args)...);
}

// Unique Ptr
template <typename T>
using Unique = std::unique_ptr<T>;

template <typename T, typename... Args>
auto MakeUnique(Args&&... args) -> decltype(std::make_unique<T>(std::forward<Args>(args)...)) {
    return std::make_unique<T>(std::forward<Args>(args)...);
}

namespace Ptr {
    /// @brief Swap function ptrs.
    /// @tparam T Type
    /// @param a Pointer 1.
    /// @param b Poitner 2.
    template <typename T>
    static inline void swap(T& a, T& b) {
        T tmp = std::move(a);
        a = std::move(b);
        b = std::move(tmp);
    }
} // End namespace Ptr

///@}

} // End namespace queso

#endif // DEFINE_INCLUDE_HPP