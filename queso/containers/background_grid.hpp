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

#pragma once

//// Project includes
#include "queso/containers/condition.hpp"
#include "queso/containers/element.hpp"
#include "queso/containers/grid_indexer.hpp"
#include "queso/containers/integration_point_concepts.hpp"
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

/// @class  BackgroundGrid
/// @author Manuel Messmer
/// @brief  Container for active elements and conditions.
///         BackgroundGrid uses GridIndexer to index the grid and to 'move' through the grid.
/// @tparam TIntegrationPoint integration-point type used by the element.
/// @tparam TBoundaryIntegrationPoint boundary integration-point type used by the element.
/// @see    grid_indexer.hpp
template<concepts::IntegrationPoint TIntegrationPoint, concepts::BoundaryIntegrationPoint TBoundaryIntegrationPoint>
class BackgroundGrid
{
public:
    ///@name Type definitions
    ///@{

    using ElementType = Element<TIntegrationPoint, TBoundaryIntegrationPoint>;
    using IntegrationPointType = TIntegrationPoint;
    using BoundaryIntegrationPointType = TBoundaryIntegrationPoint;

    using ElementContainerType = std::vector<ElementType>;
    using ElementIdMapType = std::unordered_map<IndexType, IndexType>;

    using ConditionType = Condition<ElementType>;
    using ConditionContainerType = std::vector<ConditionType>;

    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;
    using Direction = GridIndexer::Direction;

    struct NextElementResult
    {
        ElementType* pElement{};
        IndexType nextId{};
        bool isEnd{};
    };

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor
    /// @param rSettings
    explicit BackgroundGrid(const MainDictionaryType& rSettings) : mGridIndexer(rSettings)
    {}

    /// Destructor
    ~BackgroundGrid() = default;
    /// Copy Constructor
    BackgroundGrid(const BackgroundGrid& rOther) = delete;
    /// Assignment Operator
    BackgroundGrid& operator=(const BackgroundGrid& rOther) = delete;
    /// Move constructor
    BackgroundGrid(BackgroundGrid&& rOther) noexcept = default;
    /// Move assignment operator
    BackgroundGrid& operator=(BackgroundGrid&& rOther) noexcept = default;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns raw pointer to element with given id (const version).
    /// @param ElementId
    /// @return const ElementType*
    [[nodiscard]] const ElementType* pGetElement(IndexType ElementId) const
    {
        QuESo_ASSERT(mElementsLocked, "Elements must be locked (call `LockElements`) before access.");
        auto found_key = mElementIdMap.find(ElementId);
        if (found_key != mElementIdMap.end()) { return &mElements[found_key->second]; }
        return nullptr;
    }

    /// @brief Returns raw pointer to element with given id (non-const version).
    /// @param ElementId
    /// @return ElementType*
    [[nodiscard]] ElementType* pGetElement(IndexType ElementId)
    { return const_cast<ElementType*>(static_cast<const BackgroundGrid*>(this)->pGetElement(ElementId)); }

    /// @brief Returns all active elements.
    /// @return ElementContainerType&
    [[nodiscard]] ElementContainerType& GetElements() noexcept
    { return mElements; }

    /// @brief Returns all active elements.
    /// @return const ElementContainerType&
    [[nodiscard]] const ElementContainerType& GetElements() const noexcept
    { return mElements; }

    /// @brief Returns all conditions.
    /// @return const ConditionContainerType&
    [[nodiscard]] const ConditionContainerType& GetConditions() const noexcept
    { return mConditions; }

    /// @brief Adds a new condition to the background grid.
    /// @param rCondition Condition moved into the container.
    void AddCondition(ConditionType&& rCondition)
    { mConditions.emplace_back(std::move(rCondition)); }

    /// @brief Adds an element to the container.
    /// @details Requires sufficient reserved capacity to keep element addresses stable.
    /// @param rElement
    void AddElement(ElementType&& rElement)
    {
        QuESo_ASSERT(!mElementsLocked, "Cannot add elements after LockElements() was called.");
        const IndexType current_id = rElement.GetId();
        if (mElementIdMap.find(current_id) == mElementIdMap.end()) {
            const IndexType element_index = mElements.size();
            mElements.emplace_back(std::move(rElement));
            mElementIdMap.emplace(current_id, element_index);
        } else {
            QuESo_ERROR << "Element ID '" << current_id << "' already exists.\n";
        }
    }

    /// @brief Returns number of stored conditions.
    /// @return IndexType.
    [[nodiscard]] IndexType NumberOfConditions() const noexcept
    { return mConditions.size(); }

    /// @brief Returns number of stored elements.
    /// @return IndexType.
    [[nodiscard]] IndexType NumberOfActiveElements() const noexcept
    { return mElements.size(); }

    /// @brief Reserves capacity for the element container.
    /// @details Call this before inserting elements so element addresses remain stable.
    /// @param NewCapacity.
    void ReserveElements(IndexType NewCapacity)
    {
        mElements.reserve(NewCapacity);
        mElementIdMap.reserve(NewCapacity);
    }

    /// @brief Adjusts capacity of condition container.
    /// @param NewCapacity.
    void ReserveConditions(IndexType NewCapacity)
    { mConditions.reserve(NewCapacity); }

    /// @brief Locks the element container.
    /// @details After calling this function, no further elements may be added.
    ///          This guarantees stable element addresses, making pointers returned
    ///          by `pGetElement()` safe to use for the remainder of the container's lifetime.
    void LockElements() noexcept
    { mElementsLocked = true; }

    /// @brief Returns traversal information for the requested direction.
    /// @param CurrentId element id to start from.
    /// @param Dir Move direction.
    /// @return NextElementResult containing the candidate element pointer, next id, and end flag.
    [[nodiscard]] NextElementResult GetNextElement(IndexType CurrentId, Direction Dir)
    {
        switch (Dir) {
        case Direction::x_forward:
            return GetElementByIndexer(CurrentId, &GridIndexer::GetNextIndexX);
        case Direction::x_backward:
            return GetElementByIndexer(CurrentId, &GridIndexer::GetPreviousIndexX);
        case Direction::y_forward:
            return GetElementByIndexer(CurrentId, &GridIndexer::GetNextIndexY);
        case Direction::y_backward:
            return GetElementByIndexer(CurrentId, &GridIndexer::GetPreviousIndexY);
        case Direction::z_forward:
            return GetElementByIndexer(CurrentId, &GridIndexer::GetNextIndexZ);
        case Direction::z_backward:
            return GetElementByIndexer(CurrentId, &GridIndexer::GetPreviousIndexZ);
        default:
            QuESo_ASSERT(false, "There are only 6 different directions.\n");
            return {};
        }
    }

    /// @brief Returns whether the given element id is an end in the requested direction.
    /// @param CurrentId
    /// @param Dir Move direction.
    /// @return bool
    [[nodiscard]] bool IsEnd(IndexType CurrentId, Direction Dir) const
    { return mGridIndexer.IsEnd(CurrentId - 1, Dir); }

private:
    using IndexerMethodType = GridIndexer::IndexReturnType (GridIndexer::*)(IndexType) const;

    [[nodiscard]] NextElementResult GetElementByIndexer(IndexType CurrentId, IndexerMethodType Method)
    {
        const auto [next_index, index_info] = (mGridIndexer.*Method)(CurrentId - 1);
        const IndexType next_id = next_index + 1;
        bool is_end = (index_info != GridIndexer::IndexInfo::middle);
        ElementType* p_element = pGetElement(next_id);
        if (!p_element) {// Element not found. This also indicates an end.
            is_end = true;
        }
        return { p_element, next_id, is_end };
    }

    ///@}
    ///@name Private member variables
    ///@{

    GridIndexer mGridIndexer;

    ElementContainerType mElements{};
    ElementIdMapType mElementIdMap{};

    ConditionContainerType mConditions{};

    bool mElementsLocked = false;
    ///@}
};// End class BackgroundGrid
///@} // End QuESo classes

}// End namespace queso
