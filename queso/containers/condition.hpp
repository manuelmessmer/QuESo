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

//// STL includes
#include "queso/includes/settings.hpp"
#include "queso/includes/model_info.hpp"
#include "queso/containers/condition_segment.hpp"
#include "queso/containers/triangle_mesh.hpp"

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
    typedef TElementType ElementType;

    typedef ConditionSegment<ElementType> ConditionSegmentType;
    typedef Unique<ConditionSegmentType> ConditionSegmentPtrType;
    typedef std::vector<ConditionSegmentPtrType> ConditionSegmentPtrVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rConditionSettings
    /// @param rConditionInfo
    Condition( const SettingsBaseType& rConditionSettings, ModelInfoBaseType& rConditionInfo )
        :  mConditionSettings(rConditionSettings), mConditionInfo(rConditionInfo)
    {
    }

    // Destructor
    ~Condition() = default;
    // Copy constructor
    Condition(const Condition& rOther) = delete;
    // Assignement operator
    Condition& operator=(const Condition& rOther) = delete;
    /// Move constructor
    Condition(Condition&& rOther) = default;
    /// Move assignement operator
    Condition& operator=(Condition&& rOther) = default;

    /// @brief Adds new ConditionSegment to this condition. Segment is moved into container.
    /// @param pNewSegment
    void AddSegment(ConditionSegmentPtrType& pNewSegment) {
        mSegments.push_back(std::move(pNewSegment));
    }

    /// @brief Returns all stored ConditionSegments.
    /// @return const ConditionSegmentPtrVectorType&
    const ConditionSegmentPtrVectorType& GetSegments() const {
        return mSegments;
    }

    /// @brief Returns condition settings.
    /// @return const SettingsBaseType&
    const SettingsBaseType& GetSettings() const {
        return mConditionSettings;
    }

    /// @brief Returns condition settings.
    /// @return const ModelInfoBaseType&
    const ModelInfoBaseType& GetInfo() const {
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

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    DereferenceIterator<typename ConditionSegmentPtrVectorType::iterator> SegmentsBegin() {
        return dereference_iterator(mSegments.begin());
    }

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    DereferenceIterator<typename ConditionSegmentPtrVectorType::const_iterator> SegmentsBegin() const {
        return dereference_iterator(mSegments.begin());
    }

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    DereferenceIterator<typename ConditionSegmentPtrVectorType::iterator> SegmentsEnd() {
        return dereference_iterator(mSegments.end());
    }

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    DereferenceIterator<typename ConditionSegmentPtrVectorType::const_iterator> SegmentsEnd() const {
        return dereference_iterator(mSegments.end());
    }

    /// @brief Returns iterator to raw ptr. This means iterator does not point to UniquePtr<Object>,
    ///        but to Object*.
    /// @return DereferenceIterator
    RawPointerIterator<typename ConditionSegmentPtrVectorType::iterator> SegmentsBeginToPtr() {
        return raw_pointer_iterator(mSegments.begin());
    }

    /// @brief Returns iterator to raw ptr. This means iterator does not point to UniquePtr<Object>,
    ///        but to Object*.
    /// @return DereferenceIterator
    RawPointerIterator<typename ConditionSegmentPtrVectorType::const_iterator> SegmentsBeginToPtr() const {
        return raw_pointer_iterator(mSegments.begin());
    }

    /// @brief Returns iterator to raw ptr. This means iterator does not point to UniquePtr<Object>,
    ///        but to Object*.
    /// @return DereferenceIterator
    RawPointerIterator<typename ConditionSegmentPtrVectorType::iterator> SegmentsEndToPtr() {
        return raw_pointer_iterator(mSegments.end());
    }

    /// @brief Returns iterator to raw ptr. This means iterator does not point to UniquePtr<Object>,
    ///        but to Object*.
    /// @return DereferenceIterator
    RawPointerIterator<typename ConditionSegmentPtrVectorType::const_iterator> SegmentsEndToPtr() const {
        return raw_pointer_iterator(mSegments.end());
    }

private:

    ///@}
    ///@name Private Members
    ///@{

    const SettingsBaseType& mConditionSettings;
    ModelInfoBaseType& mConditionInfo;
    ConditionSegmentPtrVectorType mSegments;

    ///@}
}; // End class Condition
///@}
} // End queso namespace.

#endif // End CONDITION_INCLUDE_HPP