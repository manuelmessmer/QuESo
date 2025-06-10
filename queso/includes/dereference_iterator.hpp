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

#ifndef DEREFERENCE_ITERATOR_INCLUDE_HPP
#define DEREFERENCE_ITERATOR_INCLUDE_HPP

//// STL includes
#include <stddef.h>

namespace queso {

///@}
///@name QuESo Classes
///@{

/// @brief This iterator is intended to be used for containers like
///        std::vector<std::unique_ptr<Type>>. It allows to iterate through the container
///        without the need of dereferencing the unique_ptr.
/// @tparam BaseIterator
template <class BaseIterator>
class DereferenceIterator {
public:
    static constexpr bool is_const = std::is_const<
        std::remove_reference_t<decltype(*std::declval<BaseIterator>())>
    >::value;
    using value_type = typename BaseIterator::value_type::element_type;

    using pointer = typename std::conditional<is_const,
        const value_type*,                    // If so, use const pointer
        value_type*                           // Otherwise, use non-const pointer
    >::type;

    using reference = typename std::conditional<is_const,
        const value_type&,                      // If so, use const reference
        value_type&                             // Otherwise, use non-const reference
    >::type;

    explicit constexpr DereferenceIterator(const BaseIterator &rOther) : mIt(rOther) {}

    [[nodiscard]] constexpr reference operator*() const noexcept {
        return *(*mIt);
    }
    [[nodiscard]] constexpr pointer operator->() const noexcept {
        return mIt->get();
    }
    [[nodiscard]] constexpr reference operator[](ptrdiff_t n) const noexcept {
        return *(*(mIt + n));
    }

    constexpr DereferenceIterator& operator++() noexcept {
        ++mIt;
        return *this;
    }
    constexpr DereferenceIterator operator++(int) noexcept {
        DereferenceIterator tmp(*this);
        ++mIt;
        return tmp;
    }
    constexpr DereferenceIterator& operator+=(ptrdiff_t inc) noexcept {
        mIt += inc;
        return *this;
    }

    // Provide access to the underlying iterator for operator+
    [[nodiscard]] constexpr const BaseIterator& base() const noexcept {
        return mIt;
    }

private:
    BaseIterator mIt;
};

template <typename Iterator>
[[nodiscard]] constexpr DereferenceIterator<Iterator> dereference_iterator(Iterator t) noexcept {
    return DereferenceIterator<Iterator>(t);
}

template <class BaseIterator>
[[nodiscard]] constexpr DereferenceIterator<BaseIterator> operator+(
        const DereferenceIterator<BaseIterator>& it, ptrdiff_t n) noexcept {
    return DereferenceIterator<BaseIterator>(it.base() + n);
}

template <class BaseIterator>
[[nodiscard]] constexpr DereferenceIterator<BaseIterator> operator+(
        ptrdiff_t n, const DereferenceIterator<BaseIterator>& it) noexcept {
    return DereferenceIterator<BaseIterator>(n + it.base());
}

/// @brief Enables to use range-based for loops with DereferenceIterator.
/// @tparam TIterator
template<class TIterator>
struct DereferenceRange {
    using iterator = TIterator;

    iterator mBegin, mEnd;

    [[nodiscard]] constexpr iterator begin() const noexcept { return mBegin; }
    [[nodiscard]] constexpr iterator end() const noexcept { return mEnd; }
};

///@}

} // End namespace queso

#endif // DEREFERENCE_ITERATOR_INCLUDE_HPP