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
#include <iterator>
#include <stddef.h>

namespace queso {

///@}
///@name QuESo Classes
///@{

/// @brief This iterator is intended to be used for containers like
///        std::vector<std::unique_ptr<Type>>. It allows to iterate through the container
///        without the need of dereferencing the unique_ptr.
/// @tparam TBaseIterator (must be random access iterator).
template <class TBaseIterator>
class DereferenceIterator {
private:
    ///@}
    ///@name Helper type traits
    ///@{

    template <typename, typename = void>
    struct has_element_type : std::false_type {};

    template <typename T>
    struct has_element_type<T, std::void_t<typename T::element_type>> : std::true_type {};

    static constexpr bool it_is_const = std::is_const_v<
        std::remove_reference_t<decltype(*std::declval<TBaseIterator>())>>;

public:
    ///@}
    ///@name Static asserts
    ///@{

    static_assert( std::is_same_v<typename std::iterator_traits<TBaseIterator>::iterator_category,
                   std::random_access_iterator_tag>,
        "TBaseIterator must be a random access iterator." );

    static_assert( std::is_class_v<typename TBaseIterator::value_type>,
        "TBaseIterator::value_type must be a class (smart pointer)." );

    static_assert( has_element_type<typename TBaseIterator::value_type>::value,
        "TBaseIterator::value_type must have nested type 'element_type'." );

    ///@}
    ///@name Type definitions
    ///@{

    using iterator_category = std::random_access_iterator_tag;
    using value_type = typename TBaseIterator::value_type::element_type;
    using pointer = typename std::conditional_t<it_is_const, const value_type*, value_type*>;
    using reference = typename std::conditional_t<it_is_const, const value_type&, value_type&>;
    using difference_type = ptrdiff_t;

    ///@}
    ///@name Life cycle
    ///@{

    constexpr DereferenceIterator() noexcept = default;
    explicit constexpr DereferenceIterator( const TBaseIterator &rOther );

    template <class OtherIterator,
        std::enable_if_t<std::is_convertible_v<OtherIterator, TBaseIterator>, int> = 0>
    constexpr DereferenceIterator( const DereferenceIterator<OtherIterator>& rOther ) noexcept;

    constexpr DereferenceIterator( const DereferenceIterator& ) = default;
    constexpr DereferenceIterator( DereferenceIterator&& ) = default;
    constexpr DereferenceIterator& operator=( const DereferenceIterator& ) = default;
    constexpr DereferenceIterator& operator=( DereferenceIterator&& ) = default;

    ///@}
    ///@name Operations
    ///@{

    [[nodiscard]] constexpr reference operator[]( difference_type Inc ) const noexcept;
    [[nodiscard]] constexpr reference operator*() const noexcept;
    [[nodiscard]] constexpr pointer operator->() const noexcept;

    constexpr DereferenceIterator& operator++() noexcept;
    constexpr DereferenceIterator operator++(int) noexcept;
    constexpr DereferenceIterator& operator--() noexcept;
    constexpr DereferenceIterator operator--(int) noexcept;

    constexpr DereferenceIterator& operator+=( difference_type Inc ) noexcept;
    constexpr DereferenceIterator& operator-=( difference_type Dec ) noexcept;

    [[nodiscard]] constexpr const TBaseIterator& base() const noexcept;

private:
    ///@}
    ///@name Private members
    ///@{

    TBaseIterator mIt;

    ///@}
};

// Life cycle
template<class TBaseIterator>
constexpr DereferenceIterator<TBaseIterator>::DereferenceIterator( const TBaseIterator &rOther )
    : mIt( rOther ) {
}

template<class TBaseIterator>
template<class TOtherIterator, std::enable_if_t<std::is_convertible_v<TOtherIterator, TBaseIterator>, int>>
constexpr DereferenceIterator<TBaseIterator>::DereferenceIterator( const DereferenceIterator<TOtherIterator>& rOther ) noexcept
    : mIt( rOther.base() ) {
}

// Member operations
template<class TBaseIterator>
constexpr auto DereferenceIterator<TBaseIterator>::operator[]( difference_type Inc ) const noexcept -> reference {
    return *(*(mIt + Inc));
}

template<class TBaseIterator>
constexpr auto DereferenceIterator<TBaseIterator>::operator*() const noexcept -> reference {
    return *(*mIt);
}

template<class TBaseIterator>
constexpr auto DereferenceIterator<TBaseIterator>::operator->() const noexcept -> pointer {
    return mIt->get();
}


template<class TBaseIterator>
constexpr DereferenceIterator<TBaseIterator>& DereferenceIterator<TBaseIterator>::operator++() noexcept {
    ++mIt;
    return *this;
}

template<class TBaseIterator>
constexpr DereferenceIterator<TBaseIterator> DereferenceIterator<TBaseIterator>::operator++(int) noexcept {
    DereferenceIterator tmp(*this);
    ++mIt;
    return tmp;
}

template<class TBaseIterator>
constexpr DereferenceIterator<TBaseIterator>& DereferenceIterator<TBaseIterator>::operator--() noexcept {
    --mIt;
    return *this;
}

template<class TBaseIterator>
constexpr DereferenceIterator<TBaseIterator> DereferenceIterator<TBaseIterator>::operator--(int) noexcept {
    DereferenceIterator tmp(*this);
    --mIt;
    return tmp;
}


template<class TBaseIterator>
constexpr DereferenceIterator<TBaseIterator>&
    DereferenceIterator<TBaseIterator>::operator+=( difference_type Inc ) noexcept
{
    mIt += Inc;
    return *this;
}

template<class TBaseIterator>
constexpr DereferenceIterator<TBaseIterator>&
    DereferenceIterator<TBaseIterator>::operator-=( difference_type Dec ) noexcept
{
    mIt -= Dec;
    return *this;
}


template<class TBaseIterator>
constexpr const TBaseIterator& DereferenceIterator<TBaseIterator>::base() const noexcept {
    return mIt;
}

// Non-member operators and utility functions
template <class TBaseIterator>
[[nodiscard]] constexpr bool
    operator==(const DereferenceIterator<TBaseIterator>& rLhs, const DereferenceIterator<TBaseIterator>& rRhs) noexcept
{
    return rLhs.base() == rRhs.base();
}

template <class TBaseIterator>
[[nodiscard]] constexpr bool
    operator!=(const DereferenceIterator<TBaseIterator>& rLhs, const DereferenceIterator<TBaseIterator>& rRhs) noexcept
{
    return rLhs.base() != rRhs.base();
}

template <class TBaseIterator>
[[nodiscard]] constexpr
    bool operator<( const DereferenceIterator<TBaseIterator>& rLhs, const DereferenceIterator<TBaseIterator>& rRhs ) noexcept
{
    return rLhs.base() < rRhs.base();
}

template <class TBaseIterator>
[[nodiscard]] constexpr bool
    operator<=( const DereferenceIterator<TBaseIterator>& rLhs, const DereferenceIterator<TBaseIterator>& rRhs ) noexcept
{
    return rLhs.base() <= rRhs.base();
}

template <class TBaseIterator>
[[nodiscard]] constexpr bool
    operator>( const DereferenceIterator<TBaseIterator>& rLhs, const DereferenceIterator<TBaseIterator>& rRhs) noexcept
{
    return rLhs.base() > rRhs.base();
}

template <class TBaseIterator>
[[nodiscard]] constexpr bool
    operator>=( const DereferenceIterator<TBaseIterator>& rLhs, const DereferenceIterator<TBaseIterator>& rRhs) noexcept
{
    return rLhs.base() >= rRhs.base();
}

template <class TBaseIterator>
[[nodiscard]] constexpr DereferenceIterator<TBaseIterator>
    operator+( const DereferenceIterator<TBaseIterator>& rIt, ptrdiff_t Inc ) noexcept
{
    return DereferenceIterator<TBaseIterator>(rIt.base() + Inc);
}

template <class TBaseIterator>
[[nodiscard]] constexpr DereferenceIterator<TBaseIterator>
    operator+( ptrdiff_t Inc, const DereferenceIterator<TBaseIterator>& rIt ) noexcept
{
    return DereferenceIterator<TBaseIterator>(Inc + rIt.base());
}

template <class TBaseIterator>
[[nodiscard]] constexpr DereferenceIterator<TBaseIterator>
    operator-( const DereferenceIterator<TBaseIterator>& rIt, ptrdiff_t Dec ) noexcept
{
    return DereferenceIterator<TBaseIterator>(rIt.base() - Dec);
}

template <class TBaseIterator>
[[nodiscard]] constexpr ptrdiff_t
    operator-( const DereferenceIterator<TBaseIterator>& rLhs, const DereferenceIterator<TBaseIterator>& rRhs ) noexcept
{
    return rLhs.base() - rRhs.base();
}

/// @brief Helper function to create dereference iterator.
/// @tparam TBaseIterator
/// @param It
/// @return DereferenceIterator<TBaseIterator>
template <typename TBaseIterator>
[[nodiscard]] constexpr DereferenceIterator<TBaseIterator> dereference_iterator(TBaseIterator It) noexcept {
    return DereferenceIterator<TBaseIterator>(It);
}

/// @brief Enables to use range-based for loops with DereferenceIterator.
/// @tparam TBaseIterator
template<class TBaseIterator>
struct DereferenceRange {
    using iterator = DereferenceIterator<TBaseIterator>;

    iterator mBegin, mEnd;

    [[nodiscard]] constexpr iterator begin() const noexcept { return mBegin; }
    [[nodiscard]] constexpr iterator end() const noexcept { return mEnd; }
};

/// @brief Helper function to create dereference range.
/// @tparam TBaseIterator
/// @param It
/// @return DereferenceRange<TBaseIterator>
template <typename TBaseIterator>
[[nodiscard]] constexpr DereferenceRange<TBaseIterator>
    dereference_range(TBaseIterator Begin, TBaseIterator End) noexcept
{
    return DereferenceRange<TBaseIterator>{dereference_iterator(Begin), dereference_iterator(End)};
}

///@}

} // End namespace queso

#endif // DEREFERENCE_ITERATOR_INCLUDE_HPP