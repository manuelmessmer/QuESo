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

//// STL includes
#include <optional>
#include <span>

//// Project includes
#include "queso/containers/condition.hpp"
#include "queso/containers/element_view.hpp"
#include "queso/containers/grid_indexer.hpp"
#include "queso/containers/integration_point_concepts.hpp"
#include "queso/containers/trimmed_element.hpp"
#include "queso/containers/untrimmed_element.hpp"
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

    using UntrimmedElementType = UntrimmedElement<TIntegrationPoint, TBoundaryIntegrationPoint>;
    using TrimmedElementType = TrimmedElement<TIntegrationPoint, TBoundaryIntegrationPoint>;
    using ElementViewType = ElementView<TIntegrationPoint, TBoundaryIntegrationPoint>;
    using IntegrationPointType = TIntegrationPoint;
    using BoundaryIntegrationPointType = TBoundaryIntegrationPoint;

    using UntrimmedElementContainerType = std::vector<UntrimmedElementType>;
    using TrimmedElementContainerType = std::vector<TrimmedElementType>;

    using ConditionType = Condition<ElementViewType>;
    using ConditionContainerType = std::vector<ConditionType>;

    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;
    using Direction = GridIndexer::Direction;

    enum class ElementFilter { all, trimmed, untrimmed };

    /// Result type for view-returning traversal (always std::optional<ElementViewType>).
    struct NextElementViewResult
    {
        std::optional<ElementViewType> element{};
        IndexType next_id{};
        bool is_end{};
    };

    /// Result type for ptr-returning traversal (ElementFilter::all is not valid).
    template<ElementFilter TFilter>
    struct NextElementResult
    {
        static_assert(
            TFilter != ElementFilter::all,
            "ElementFilter::all is not valid for NextElementResult — use NextElementViewResult instead."
        );
        using ElementPtrType =
            std::conditional_t<TFilter == ElementFilter::untrimmed, UntrimmedElementType*, TrimmedElementType*>;
        ElementPtrType p_element{};
        IndexType next_id{};
        bool is_end{};
    };

    /// @brief Describes the location of an element within the background-grid containers.
    struct ElementStorageLocation
    {
        bool is_trimmed{};
        IndexType index{};
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

    /// @brief Returns an ElementView for the element with the given id.
    /// @details TFilter controls which elements are considered:
    ///          all       → searches both containers.
    ///          trimmed   → only returns a value if the element is trimmed.
    ///          untrimmed → only returns a value if the element is untrimmed.
    /// @param ElementId
    /// @return std::optional<ElementViewType>
    template<ElementFilter TFilter = ElementFilter::all>
    [[nodiscard]] std::optional<ElementViewType> GetElementView(IndexType ElementId) const
    {
        QuESo_ASSERT(mElementsLocked, "Elements must be locked (call `LockElements`) before access.");
        auto found_key = mElementIdMap.find(ElementId);
        if (found_key == mElementIdMap.end()) { return std::nullopt; }
        const auto& rLocation = found_key->second;
        if constexpr (TFilter == ElementFilter::all) {
            if (rLocation.is_trimmed) return mTrimmedElements[rLocation.index].View();
            return mUntrimmedElements[rLocation.index].View();
        } else if constexpr (TFilter == ElementFilter::untrimmed) {
            if (!rLocation.is_trimmed) return mUntrimmedElements[rLocation.index].View();
            return std::nullopt;
        } else {
            if (rLocation.is_trimmed) return mTrimmedElements[rLocation.index].View();
            return std::nullopt;
        }
    }

    /// @brief Returns a mutable pointer to the element with the given id.
    /// @details ElementFilter::all is not valid — use GetElementView() instead.
    ///          Returns nullptr if the element does not exist or does not match TFilter.
    /// @param ElementId
    template<ElementFilter TFilter>
    [[nodiscard]] decltype(auto) pGetElement(IndexType ElementId)
    {
        static_assert(
            TFilter != ElementFilter::all,
            "ElementFilter::all is not valid for pGetElement — use GetElementView() instead."
        );
        QuESo_ASSERT(mElementsLocked, "Elements must be locked (call `LockElements`) before access.");

        auto found_key = mElementIdMap.find(ElementId);
        if constexpr (TFilter == ElementFilter::untrimmed) {
            if (found_key != mElementIdMap.end() && !found_key->second.is_trimmed)
                return &mUntrimmedElements[found_key->second.index];
            return static_cast<UntrimmedElementType*>(nullptr);
        } else {
            if (found_key != mElementIdMap.end() && found_key->second.is_trimmed)
                return &mTrimmedElements[found_key->second.index];
            return static_cast<TrimmedElementType*>(nullptr);
        }
    }

    /// @brief Returns a lazy range of ElementView over all matching elements.
    /// @details TFilter controls which elements are included:
    ///          all       → both trimmed and untrimmed elements.
    ///          untrimmed → only untrimmed elements.
    ///          trimmed   → only trimmed elements.
    /// @return lazy range of ElementViewType
    template<ElementFilter TFilter = ElementFilter::all>
    [[nodiscard]] auto GetElementViews() const
    {
        if constexpr (TFilter == ElementFilter::all) {
            const auto num_untrimmed = mUntrimmedElements.size();
            const auto num_total = num_untrimmed + mTrimmedElements.size();
            return std::views::iota(SizeType{ 0 }, num_total)
                   | std::views::transform([this, num_untrimmed](SizeType i) {
                         if (i < num_untrimmed) return mUntrimmedElements[i].View();
                         return mTrimmedElements[i - num_untrimmed].View();
                     });
        } else if constexpr (TFilter == ElementFilter::untrimmed) {
            return std::views::iota(SizeType{ 0 }, mUntrimmedElements.size())
                   | std::views::transform([this](SizeType i) { return mUntrimmedElements[i].View(); });
        } else {
            return std::views::iota(SizeType{ 0 }, mTrimmedElements.size())
                   | std::views::transform([this](SizeType i) { return mTrimmedElements[i].View(); });
        }
    }

    /// @brief Returns a span over the raw element container matching TFilter.
    /// @details ElementFilter::all is not valid (the two containers have different element types).
    ///          Returns std::span<const UntrimmedElementType> or std::span<const TrimmedElementType>.
    template<ElementFilter TFilter>
    [[nodiscard]] auto GetElements() const noexcept
    {
        static_assert(
            TFilter != ElementFilter::all,
            "ElementFilter::all is not valid for GetElements — use GetElementViews() instead."
        );
        if constexpr (TFilter == ElementFilter::untrimmed)
            return std::span<const UntrimmedElementType>(mUntrimmedElements);
        else
            return std::span<const TrimmedElementType>(mTrimmedElements);
    }

    /// @brief Returns a mutable span over the raw element container matching TFilter.
    /// @details ElementFilter::all is not valid (the two containers have different element types).
    ///          Returns std::span<UntrimmedElementType> or std::span<TrimmedElementType>.
    template<ElementFilter TFilter>
    [[nodiscard]] auto GetElements() noexcept
    {
        static_assert(
            TFilter != ElementFilter::all,
            "ElementFilter::all is not valid for GetElements — use GetElementViews() instead."
        );
        if constexpr (TFilter == ElementFilter::untrimmed)
            return std::span<UntrimmedElementType>(mUntrimmedElements);
        else
            return std::span<TrimmedElementType>(mTrimmedElements);
    }

    /// @brief Returns all conditions.
    /// @return const ConditionContainerType&
    [[nodiscard]] const ConditionContainerType& GetConditions() const noexcept
    { return mConditions; }

    /// @brief Adds a new condition to the background grid.
    /// @param rCondition Condition moved into the container.
    void AddCondition(ConditionType&& rCondition)
    { mConditions.emplace_back(std::move(rCondition)); }

    /// @brief Returns number of stored conditions.
    /// @return IndexType.
    [[nodiscard]] IndexType NumberOfConditions() const noexcept
    { return mConditions.size(); }

    /// @brief Returns number of stored elements (trimmed + untrimmed).
    /// @return IndexType.
    [[nodiscard]] IndexType NumberOfActiveElements() const noexcept
    { return mUntrimmedElements.size() + mTrimmedElements.size(); }

    /// @brief Returns number of stored untrimmed elements.
    /// @return IndexType.
    [[nodiscard]] IndexType NumberOfUntrimmedElements() const noexcept
    { return mUntrimmedElements.size(); }

    /// @brief Returns number of stored trimmed elements.
    /// @return IndexType.
    [[nodiscard]] IndexType NumberOfTrimmedElements() const noexcept
    { return mTrimmedElements.size(); }

    /// @brief Reserves capacity for the element container.
    /// @details Call this before inserting elements so element addresses remain stable.
    /// @param NewCapacity.
    void ReserveElements(IndexType NewCapacity)
    {
        mElementIdMap.reserve(NewCapacity);
        mUntrimmedElements.reserve(NewCapacity);
        mTrimmedElements.reserve(NewCapacity);
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

    /// @brief Traverses the grid in Dir and returns an ElementView of the next element.
    /// @details TFilter controls which elements are considered (all / trimmed / untrimmed).
    ///          Returns NextElementViewResult with std::optional<ElementViewType>.
    /// @param CurrentId element id to start from.
    /// @param Dir Move direction.
    template<ElementFilter TFilter = ElementFilter::all>
    [[nodiscard]] NextElementViewResult GetNextElementView(IndexType CurrentId, Direction Dir) const
    {
        switch (Dir) {
        case Direction::x_forward:
            return GetElementViewByIndexer<TFilter>(CurrentId, &GridIndexer::GetNextIndexX);
        case Direction::x_backward:
            return GetElementViewByIndexer<TFilter>(CurrentId, &GridIndexer::GetPreviousIndexX);
        case Direction::y_forward:
            return GetElementViewByIndexer<TFilter>(CurrentId, &GridIndexer::GetNextIndexY);
        case Direction::y_backward:
            return GetElementViewByIndexer<TFilter>(CurrentId, &GridIndexer::GetPreviousIndexY);
        case Direction::z_forward:
            return GetElementViewByIndexer<TFilter>(CurrentId, &GridIndexer::GetNextIndexZ);
        case Direction::z_backward:
            return GetElementViewByIndexer<TFilter>(CurrentId, &GridIndexer::GetPreviousIndexZ);
        default:
            QuESo_ASSERT(false, "There are only 6 different directions.\n");
            return {};
        }
    }

    /// @brief Traverses the grid in Dir and returns a mutable pointer to the next element.
    /// @details ElementFilter::all is not valid — use GetNextElementView() instead.
    ///          Returns NextElementResult<TFilter> with a raw pointer (nullptr if not found).
    /// @param CurrentId element id to start from.
    /// @param Dir Move direction.
    template<ElementFilter TFilter>
    [[nodiscard]] NextElementResult<TFilter> GetNextElement(IndexType CurrentId, Direction Dir)
    {
        static_assert(
            TFilter != ElementFilter::all,
            "ElementFilter::all is not valid for GetNextElement — use GetNextElementView() instead."
        );
        switch (Dir) {
        case Direction::x_forward:
            return GetElementByIndexer<TFilter>(CurrentId, &GridIndexer::GetNextIndexX);
        case Direction::x_backward:
            return GetElementByIndexer<TFilter>(CurrentId, &GridIndexer::GetPreviousIndexX);
        case Direction::y_forward:
            return GetElementByIndexer<TFilter>(CurrentId, &GridIndexer::GetNextIndexY);
        case Direction::y_backward:
            return GetElementByIndexer<TFilter>(CurrentId, &GridIndexer::GetPreviousIndexY);
        case Direction::z_forward:
            return GetElementByIndexer<TFilter>(CurrentId, &GridIndexer::GetNextIndexZ);
        case Direction::z_backward:
            return GetElementByIndexer<TFilter>(CurrentId, &GridIndexer::GetPreviousIndexZ);
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

    /// @brief Calls rBuilder.Build(), then inserts the element into trimmed storage.
    /// @details Dispatches to this overload when rBuilder.Build() returns std::optional<T>,
    ///          i.e., when the builder may reject the element (no domain or no IPs).
    ///          The critical section protects concurrent writes from the OMP parallel loop.
    /// @tparam TBuilderType  Must provide Build(IndexType, const ElementBounds&)
    ///                        returning std::optional<TrimmedElement>.
    /// @return true if the element was built and inserted successfully.
    template<typename TBuilderType>
        requires requires(TBuilderType& b, IndexType id, const ElementBounds& bounds) { b.Build(id, bounds).value(); }
    bool MakeElement(TBuilderType& rBuilder, IndexType Id, const ElementBounds& rBounds)
    {
        auto result = rBuilder.Build(Id, rBounds);
        if (!result) return false;
        QuESo_ASSERT(!mElementsLocked, "Cannot add elements after LockElements() was called.");
#pragma omp critical
        {
            QuESo_ERROR_IF(mElementIdMap.contains(Id)) << "Element ID '" << Id << "' already exists.\n";
            const IndexType element_index = mTrimmedElements.size();
            mTrimmedElements.emplace_back(std::move(*result));
            mElementIdMap.emplace(Id, ElementStorageLocation{ true, element_index });
        }
        return true;
    }

    /// @brief Calls rBuilder.Build(), then inserts the element into untrimmed storage.
    /// @details Dispatches to this overload when rBuilder.Build() returns T directly
    ///          (non-optional), i.e., the builder always produces a valid element.
    ///          The critical section protects concurrent writes from the OMP parallel loop.
    /// @tparam TBuilderType  Must provide Build(IndexType, const ElementBounds&)
    ///                        returning UntrimmedElement.
    template<typename TBuilderType>
        requires(!requires(TBuilderType& b, IndexType id, const ElementBounds& bounds) { b.Build(id, bounds).value(); })
    void MakeElement(TBuilderType& rBuilder, IndexType Id, const ElementBounds& rBounds)
    {
        QuESo_ASSERT(!mElementsLocked, "Cannot add elements after LockElements() was called.");
        auto element = rBuilder.Build(Id, rBounds);
#pragma omp critical
        {
            QuESo_ERROR_IF(mElementIdMap.contains(Id)) << "Element ID '" << Id << "' already exists.\n";
            const IndexType element_index = mUntrimmedElements.size();
            mUntrimmedElements.emplace_back(std::move(element));
            mElementIdMap.emplace(Id, ElementStorageLocation{ false, element_index });
        }
    }

private:
    using IndexerMethodType = GridIndexer::IndexReturnType (GridIndexer::*)(IndexType) const;

    template<ElementFilter TFilter = ElementFilter::all>
    [[nodiscard]] NextElementViewResult GetElementViewByIndexer(IndexType CurrentId, IndexerMethodType Method) const
    {
        const auto [next_index, index_info] = (mGridIndexer.*Method)(CurrentId - 1);
        const IndexType next_id = next_index + 1;
        bool is_end = (index_info != GridIndexer::IndexInfo::middle);
        auto p_element = GetElementView<TFilter>(next_id);
        if (!p_element) is_end = true;
        return { p_element, next_id, is_end };
    }

    template<ElementFilter TFilter>
    [[nodiscard]] NextElementResult<TFilter> GetElementByIndexer(IndexType CurrentId, IndexerMethodType Method)
    {
        static_assert(
            TFilter != ElementFilter::all,
            "ElementFilter::all is not valid for GetElementByIndexer — use GetElementViewByIndexer() instead."
        );
        const auto [next_index, index_info] = (mGridIndexer.*Method)(CurrentId - 1);
        const IndexType next_id = next_index + 1;
        bool is_end = (index_info != GridIndexer::IndexInfo::middle);
        auto p_element = pGetElement<TFilter>(next_id);
        if (!p_element) is_end = true;
        return { p_element, next_id, is_end };
    }

    ///@}
    ///@name Private member variables
    ///@{

    GridIndexer mGridIndexer;
    UntrimmedElementContainerType mUntrimmedElements{};
    TrimmedElementContainerType mTrimmedElements{};
    std::unordered_map<IndexType, ElementStorageLocation> mElementIdMap{};
    ConditionContainerType mConditions{};
    bool mElementsLocked = false;

    ///@}
};// End class BackgroundGrid
///@} // End QuESo classes

}// End namespace queso
