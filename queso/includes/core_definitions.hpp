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

#ifndef CORE_DEFINITIONS_INCLUDE_HPP
#define CORE_DEFINITIONS_INCLUDE_HPP

//// STL includes
#include <ostream>
#include <limits>
#include <memory>
#include <array>

namespace queso {

///@name QuESo GLOBAL TYPE DEFINITIONS
///@{

typedef std::size_t  SizeType;
typedef std::size_t  IndexType;

typedef std::array<double,3> PointType;
typedef std::array<double,3> Vector3d;
typedef std::array<IndexType,3>  Vector3i;
typedef std::pair<PointType, PointType> BoundingBoxType;
typedef std::pair<Vector3i, Vector3i> PartitionBoxType;

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
    const PointType delta{rUpperBound[0] - rLowerBound[0],
                          rUpperBound[1] - rLowerBound[1],
                          rUpperBound[2] - rLowerBound[2]};
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
enum class IntegrationMethod {gauss, gauss_reduced_1, gauss_reduced_2, ggq_optimal, ggq_reduced_1, ggq_reduced_2};
typedef IntegrationMethod IntegrationMethodType;
inline std::ostream& operator<<(std::ostream& rOs, IntegrationMethodType Enum) {
    switch(Enum) {
        case IntegrationMethod::gauss:
            return (rOs << "Gauss");
        case IntegrationMethod::gauss_reduced_1:
            return (rOs << "Gauss_Reduced1");
        case IntegrationMethod::gauss_reduced_2:
            return (rOs << "Gauss_Reduced2");
        case IntegrationMethod::ggq_optimal:
            return (rOs << "GGQ_Optimal");
        case IntegrationMethod::ggq_reduced_1:
            return (rOs << "GGQ_Reduced1");
        case IntegrationMethod::ggq_reduced_2:
            return (rOs << "GGQ_Reduced2");
        default:
            return rOs;
    }
}

enum class IntersectionState {inside, outside, trimmed};
typedef IntersectionState IntersectionStateType;
inline std::ostream& operator<<(std::ostream& rOs, IntersectionStateType Enum) {
    switch(Enum) {
        case IntersectionState::inside:
            return (rOs << "Inside");
        case IntersectionState::outside:
            return (rOs << "Outside");
        case IntersectionState::trimmed:
            return (rOs << "Trimmed");
        default:
            return rOs;
    }
}

enum class GridType {b_spline_grid, hexahedral_fe_grid};
typedef GridType GridTypeType;
inline std::ostream& operator<<(std::ostream& rOs, GridTypeType Enum) {
    switch(Enum) {
        case GridType::b_spline_grid:
            return (rOs << "b_spline_grid");
        case GridType::hexahedral_fe_grid:
            return (rOs << "hexahedral_fe_grid");
        default:
            return rOs;
    }
}

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

///@}
///@name QuESo ITERATOR DEFINITIONS
///@{

// DereferenceIterator
template <class BaseIterator> class DereferenceIterator : public BaseIterator {
public:
    using value_type = typename BaseIterator::value_type::element_type;
    using pointer = value_type *;
    using reference = value_type &;

    DereferenceIterator(const BaseIterator &other) : BaseIterator(other) {}

    reference operator*() const { return *(this->BaseIterator::operator*()); }
    pointer operator->() const { return this->BaseIterator::operator*().get(); }
    reference operator[](size_t n) const {
        return *(this->BaseIterator::operator[](n));
    }
};

template <typename Iterator> DereferenceIterator<Iterator> dereference_iterator(Iterator t) {
    return DereferenceIterator<Iterator>(t);
}

// RawPointerIterator
template <class BaseIterator> class RawPointerIterator : public BaseIterator {
public:
    using value_type = typename BaseIterator::value_type::element_type;
    using pointer = value_type *;
    using reference = value_type &;

    RawPointerIterator(const BaseIterator &other) : BaseIterator(other) {}

    pointer operator*() const { return (this->BaseIterator::operator*()).get(); }
    pointer operator->() const { return this->BaseIterator::operator*().get(); }
    pointer operator[](size_t n) const {
        return (this->BaseIterator::operator[](n)).get();
    }
};

template <typename Iterator> RawPointerIterator<Iterator> raw_pointer_iterator(Iterator t) {
    return RawPointerIterator<Iterator>(t);
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


template<typename type>
inline std::ostream& operator<<(std::ostream& rOStream, const std::array<type, 3>& rThis)  {
    rOStream << '(' << rThis[0] << ", " << rThis[1] << ", " << rThis[2] << ')';
    return rOStream;
}
///@}

} // End namespace queso

#endif // CORE_DEFINITIONS_INCLUDE_HPP