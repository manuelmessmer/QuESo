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
        const auto& r_location = found_key->second;
        if constexpr (TFilter == ElementFilter::all) {
            if (r_location.is_trimmed) return mTrimmedElements[r_location.index].View();
            return mUntrimmedElements[r_location.index].View();
        } else if constexpr (TFilter == ElementFilter::untrimmed) {
            if (!r_location.is_trimmed) return mUntrimmedElements[r_location.index].View();
            return std::nullopt;
        } else {
            if (r_location.is_trimmed) return mTrimmedElements[r_location.index].View();
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

    /// @brief Reserves capacity for the element containers.
    /// @param NewCapacity.
    void ReserveElements(IndexType NumUntrimmedElements, IndexType NumTrimmedElements)
    {
        mElementIdMap.reserve(NumUntrimmedElements + NumTrimmedElements);

        mUntrimmedElements.reserve(NumUntrimmedElements);
        mTrimmedElements.reserve(NumTrimmedElements);
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

    /// @brief Traverses the grid in the requested direction and returns a view of the next element.
    /// @details TFilter controls which elements are considered:
    ///          `all`, `trimmed`, or `untrimmed`.
    /// @param CurrentId Element id to start from.
    /// @param Dir Move direction.
    /// @return NextElementViewResult containing the optional next element, the next grid id, and the end flag.
    template<ElementFilter TFilter = ElementFilter::all>
    [[nodiscard]] NextElementViewResult GetNextElementView(IndexType CurrentId, Direction Dir) const
    {
        return GetNextElementByDirection(CurrentId, Dir, [this](IndexType NextId) {
            return GetElementView<TFilter>(NextId);
        });
    }

    /// @brief Traverses the grid in the requested direction and returns a mutable pointer to the next element.
    /// @details `ElementFilter::all` is not valid here; use GetNextElementView() instead.
    /// @param CurrentId Element id to start from.
    /// @param Dir Move direction.
    /// @return NextElementResult<TFilter> containing the pointer to the next element, the next grid id, and the end
    ///         flag.
    template<ElementFilter TFilter>
    [[nodiscard]] NextElementResult<TFilter> GetNextElement(IndexType CurrentId, Direction Dir)
    {
        static_assert(
            TFilter != ElementFilter::all,
            "ElementFilter::all is not valid for GetNextElement — use GetNextElementView() instead."
        );
        return GetNextElementByDirection(CurrentId, Dir, [this](IndexType NextId) {
            return pGetElement<TFilter>(NextId);
        });
    }

    /// @brief Returns whether the given element id is an end in the requested direction.
    /// @param CurrentId
    /// @param Dir Move direction.
    /// @return bool
    [[nodiscard]] bool IsEnd(IndexType CurrentId, Direction Dir) const
    { return mGridIndexer.IsEnd(CurrentId - 1, Dir); }

    /// @brief Calls rBuilder.Build() and inserts the element into the corresponding storage.
    /// @details The target storage is selected from `TBuilderType::Builds` (e.g., ElementFilter::trimmed).
    ///          If Build() returns an empty result, no element is inserted.
    ///          The critical section protects concurrent writes from the OpenMP loop.
    /// @return true if an element was built and inserted, false otherwise.
    template<typename TBuilderType>
    bool MakeElement(TBuilderType& rBuilder, IndexType Id, const ElementBounds& rBounds)
    {
        QuESo_ASSERT(!mElementsLocked, "Cannot add elements after LockElements() was called.");
        auto result = rBuilder.Build(Id, rBounds);
        if (!result) return false;
#pragma omp critical
        {
            QuESo_ERROR_IF(mElementIdMap.contains(Id)) << "Element ID '" << Id << "' already exists.\n";
            const auto location = [&]() {
                if constexpr (TBuilderType::Builds == ElementFilter::trimmed) {
                    const IndexType index = mTrimmedElements.size();
                    mTrimmedElements.emplace_back(std::move(*result));
                    return ElementStorageLocation{ true, index };
                } else {
                    static_assert(
                        TBuilderType::Builds == ElementFilter::untrimmed,
                        "The builder must either support trimmed or untrimmed element construction."
                    );

                    const IndexType index = mUntrimmedElements.size();
                    mUntrimmedElements.emplace_back(std::move(*result));
                    return ElementStorageLocation{ false, index };
                }
            }();
            mElementIdMap.emplace(Id, location);
        }
        return true;
    }


private:
    using IndexerMethodType = GridIndexer::IndexReturnType (GridIndexer::*)(IndexType) const;

    /// @brief Traverses the grid in the requested direction and dispatches to the corresponding GridIndexer step.
    /// @details The accessor decides whether the result is returned as an ElementView or as a mutable raw pointer.
    /// @param CurrentId Element id to start from.
    /// @param Dir Move direction.
    /// @param rAccessor Callable taking the next element id and returning either `std::optional<ElementViewType>`,
    ///                  `UntrimmedElementType*`, or `TrimmedElementType*`.
    template<typename TAccessor>
    [[nodiscard]] auto GetNextElementByDirection(IndexType CurrentId, Direction Dir, TAccessor&& rAccessor) const
    {
        switch (Dir) {
        case Direction::x_forward:
            return GetNextElementByIndexer(CurrentId, &GridIndexer::GetNextIndexX, std::forward<TAccessor>(rAccessor));
        case Direction::x_backward:
            return GetNextElementByIndexer(
                CurrentId, &GridIndexer::GetPreviousIndexX, std::forward<TAccessor>(rAccessor)
            );
        case Direction::y_forward:
            return GetNextElementByIndexer(CurrentId, &GridIndexer::GetNextIndexY, std::forward<TAccessor>(rAccessor));
        case Direction::y_backward:
            return GetNextElementByIndexer(
                CurrentId, &GridIndexer::GetPreviousIndexY, std::forward<TAccessor>(rAccessor)
            );
        case Direction::z_forward:
            return GetNextElementByIndexer(CurrentId, &GridIndexer::GetNextIndexZ, std::forward<TAccessor>(rAccessor));
        case Direction::z_backward:
            return GetNextElementByIndexer(
                CurrentId, &GridIndexer::GetPreviousIndexZ, std::forward<TAccessor>(rAccessor)
            );
        }
        QuESo_ERROR << "There are only 6 different directions.\n";
    }

    /// @brief Applies a GridIndexer step, fetches the corresponding element, and packages the traversal result.
    /// @details The result type is deduced from `rAccessor`:
    ///          `std::optional<ElementViewType>` yields NextElementViewResult,
    ///          `UntrimmedElementType*` yields NextElementResult<ElementFilter::untrimmed>, and
    ///          `TrimmedElementType*` yields NextElementResult<ElementFilter::trimmed>.
    /// @param CurrentId Element id to start from.
    /// @param Method GridIndexer stepping method.
    /// @param rAccessor Callable taking the computed next id and returning the matching element access type.
    template<typename TAccessor>
    [[nodiscard]] auto
        GetNextElementByIndexer(IndexType CurrentId, IndexerMethodType Method, TAccessor&& rAccessor) const
    {
        const auto [next_index, index_info] = (mGridIndexer.*Method)(CurrentId - 1);
        const IndexType next_id = next_index + 1;
        bool is_end = (index_info != GridIndexer::IndexInfo::middle);
        auto element = rAccessor(next_id);
        if (!element) { is_end = true; }

        using AccessResultType = decltype(element);
        if constexpr (std::same_as<AccessResultType, std::optional<ElementViewType>>) {
            return NextElementViewResult{ std::move(element), next_id, is_end };
        } else if constexpr (std::same_as<AccessResultType, UntrimmedElementType*>) {
            return NextElementResult<ElementFilter::untrimmed>{ element, next_id, is_end };
        } else {
            static_assert(std::same_as<AccessResultType, TrimmedElementType*>);
            return NextElementResult<ElementFilter::trimmed>{ element, next_id, is_end };
        }
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
