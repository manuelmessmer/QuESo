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
#include <functional>

//// Project includes
#include "queso/containers/condition_segment.hpp"
#include "queso/containers/dictionary.hpp"

namespace queso {

///@name QuESo Classes
///@{

/// @class  Condition
/// @author Manuel Messmer
/// @brief  Container for condition settings, condition info, and the corresponding condition segments.
///         Each segment is clipped to the element boundaries in the background grid and holds the
///         respective section of the triangle mesh.
/// @see    containers/condition_segment.h
template<typename TElementType>
class Condition
{
public:
    ///@name Type definitions
    ///@{
    using ElementType = TElementType;

    using ConditionSegmentType = ConditionSegment<ElementType>;
    using ConditionSegmentContainerType = std::vector<ConditionSegmentType>;
    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor
    /// @param rConditionSettings
    /// @param rConditionInfo
    Condition(const MainDictionaryType& rConditionSettings, MainDictionaryType& rConditionInfo)
        : mConditionSettings(rConditionSettings), mConditionInfo(rConditionInfo)
    {}

    /// Destructor
    ~Condition() = default;
    /// Copy constructor
    Condition(const Condition& rOther) = delete;
    /// Assignment operator
    Condition& operator=(const Condition& rOther) = delete;
    /// Move constructor
    Condition(Condition&& rOther) noexcept = default;
    /// Move assignment operator
    Condition& operator=(Condition&& rOther) noexcept = default;

    /// @brief Adds new ConditionSegment to this condition. Segment is moved into container.
    /// @param rNewSegment
    void AddSegment(ConditionSegmentType&& rNewSegment)
    { mSegments.emplace_back(std::move(rNewSegment)); }

    /// @brief Returns all stored ConditionSegments.
    /// @return const ConditionSegmentContainerType&
    [[nodiscard]] const ConditionSegmentContainerType& GetSegments() const noexcept
    { return mSegments; }

    /// @brief Returns condition settings.
    /// @return const MainDictionaryType&
    [[nodiscard]] const MainDictionaryType& GetSettings() const noexcept
    { return mConditionSettings.get(); }

    /// @brief Returns condition info.
    /// @return const MainDictionaryType&
    [[nodiscard]] const MainDictionaryType& GetInfo() const noexcept
    { return mConditionInfo.get(); }

    /// @brief Returns the number of segments.
    /// @return IndexType
    [[nodiscard]] IndexType NumberOfSegments() const noexcept
    { return mSegments.size(); }

private:
    ///@}
    ///@name Private Members
    ///@{

    std::reference_wrapper<const MainDictionaryType> mConditionSettings;
    std::reference_wrapper<MainDictionaryType> mConditionInfo;
    ConditionSegmentContainerType mSegments{};

    ///@}
};// End class Condition
///@}
}// namespace queso
