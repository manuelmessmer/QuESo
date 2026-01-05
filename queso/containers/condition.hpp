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

#ifndef CONDITION_INCLUDE_HPP
#define CONDITION_INCLUDE_HPP

//// Project includes
#include "queso/containers/condition_segment.hpp"
#include "queso/containers/dictionary.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  Condition
 * @author Manuel Messmer
 * @brief  Interface for conditions. Stores condition settings, some condition info, and a list of condition segments.
 *         Each segment is clipped to the element boundaries in the background grid and holds the respective
 *         section of the triangle mesh.
 * @see    containers/condition_segment.h
**/
template<typename TElementType>
class Condition {
public:

    ///@name Type Definitions
    ///@{
    using ElementType = TElementType;

    using ConditionSegmentType = ConditionSegment<ElementType>;
    using ConditionSegmentPtrType = Unique<ConditionSegmentType>;
    using ConditionSegmentContainerType = std::vector<ConditionSegmentPtrType>;
    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rConditionSettings
    /// @param rConditionInfo
    Condition( const MainDictionaryType& rConditionSettings, MainDictionaryType& rConditionInfo )
        :  mConditionSettings(rConditionSettings), mConditionInfo(rConditionInfo)
    {
    }

    /// Destructor
    ~Condition() = default;
    /// Copy constructor
    Condition(const Condition& rOther) = delete;
    /// Assignement operator
    Condition& operator=(const Condition& rOther) = delete;
    /// Move constructor
    Condition(Condition&& rOther) noexcept = default;
    /// Move assignement operator
    Condition& operator=(Condition&& rOther) noexcept = default;

    /// @brief Adds new ConditionSegment to this condition. Segment is moved into container.
    /// @param pNewSegment
    void AddSegment(ConditionSegmentPtrType& pNewSegment) {
        mSegments.push_back(std::move(pNewSegment));
    }

    /// @brief Returns all stored ConditionSegments.
    /// @return const ConditionSegmentContainerType&
    const ConditionSegmentContainerType& GetSegments() const {
        return mSegments;
    }

    /// @brief Returns condition settings.
    /// @return const SettingsBaseType&
    const MainDictionaryType& GetSettings() const {
        return mConditionSettings;
    }

    /// @brief Returns condition settings.
    /// @return const ModelInfoBaseType&
    const MainDictionaryType& GetInfo() const {
        return mConditionInfo;
    }

    /// @brief Returns the number of segments.
    /// @return IndexType
    IndexType NumberOfSegments() const {
        return mSegments.size();
    }

    ///@}
    ///@name Iterators
    ///@{

    /// @brief Returns dereference iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    auto SegmentsBegin() -> DereferenceIterator<typename ConditionSegmentContainerType::iterator> {
        return dereference_iterator(mSegments.begin());
    }

    /// @brief Returns dereference iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    auto SegmentsBegin() const -> DereferenceIterator<typename ConditionSegmentContainerType::const_iterator> {
        return dereference_iterator(mSegments.begin());
    }

    /// @brief Returns dereference iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    auto SegmentsEnd() -> DereferenceIterator<typename ConditionSegmentContainerType::iterator> {
        return dereference_iterator(mSegments.end());
    }

    /// @brief Returns dereference iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    auto SegmentsEnd() const -> DereferenceIterator<typename ConditionSegmentContainerType::const_iterator> {
        return dereference_iterator(mSegments.end());
    }

    /// @brief Returns dereference range for the segments container (non-const version).
    /// @return DereferenceRange.
    auto Segments() -> DereferenceRange<typename ConditionSegmentContainerType::iterator> {
        return dereference_range( mSegments.begin(), mSegments.end() );
    }

    /// @brief Returns dereference range for the segments container (const version).
    /// @return DereferenceRange.
    auto Segments() const -> DereferenceRange<typename ConditionSegmentContainerType::const_iterator> {
        return dereference_range( mSegments.begin(), mSegments.end() );
    }

private:

    ///@}
    ///@name Private Members
    ///@{

    const MainDictionaryType& mConditionSettings;
    MainDictionaryType& mConditionInfo;
    ConditionSegmentContainerType mSegments;

    ///@}
}; // End class Condition
///@}
} // End queso namespace.

#endif // End CONDITION_INCLUDE_HPP
