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
#include <atomic>
#include <functional>
#include <optional>

//// Project includes
#include "queso/containers/trimmed_element.hpp"
#include "queso/containers/untrimmed_element.hpp"
#include "queso/containers/background_grid.hpp"
#include "queso/embedding/brep_operator.h"
#include "queso/includes/register_keys.hpp"
#include "queso/includes/timer.hpp"
#include "queso/quadrature/single_element.hpp"
#include "queso/quadrature/trimmed_element.hpp"

namespace queso {

///@name QuESo Classes
///@{

/// @class  TrimmedElementBuilder
/// @brief Encapsulates construction policy for trimmed elements.
/// TODO: Add builder concept.
template<
    concepts::IntegrationPoint TIntegrationPointType,
    concepts::BoundaryIntegrationPoint TBoundaryIntegrationPointType>
class TrimmedElementBuilder
{
public:
    ///@name Type definitions
    ///@{
    using IntegrationPointType = TIntegrationPointType;
    using BoundaryIntegrationPointType = TBoundaryIntegrationPointType;
    using ElementType = TrimmedElement<IntegrationPointType, BoundaryIntegrationPointType>;
    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;
	using BackgroundGridType = BackgroundGrid<IntegrationPointType, BoundaryIntegrationPointType>; 
	using ElementFilterType = BackgroundGridType::ElementFilter;
	static constexpr ElementFilterType Builds = ElementFilterType::trimmed;
    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructs a TrimmedElementBuilder from the main settings dictionary.
    /// @param rSettings  The full settings dictionary. Reads trimmed quadrature,
    ///                   background-grid and general-settings sub-dictionaries.
    /// @param rBrepOperator  BRep operator used for domain extraction. Stored by reference.
    TrimmedElementBuilder(const MainDictionaryType& rSettings, const BRepOperator& rBrepOperator)
        : mBrepOperator(std::cref(rBrepOperator)),
          mMinVolRatio(
              std::max<double>(
                  rSettings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<double>(
                      TrimmedQuadratureRuleSettings::min_element_volume_ratio
                  ),
                  1e-10
              )
          ),
          mMinNumBoundaryTriangles(
              rSettings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<IndexType>(
                  TrimmedQuadratureRuleSettings::min_num_boundary_triangles
              )
          ),
          mNeglectIfStlIsFlawed(
              rSettings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<bool>(
                  TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed
              )
          ),
          mPolynomialOrder(
              rSettings[MainSettings::background_grid_settings].GetRequiredValue<Vector3i>(
                  BackgroundGridSettings::polynomial_order
              )
          ),
          mMomentFittingResidual(
              rSettings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<double>(
                  TrimmedQuadratureRuleSettings::moment_fitting_residual
              )
          ),
          mEchoLevel(rSettings[MainSettings::general_settings].GetRequiredValue<IndexType>(GeneralSettings::echo_level))
    {}

    ///@}
    ///@name Operations
    ///@{

    /// @brief Builds a trimmed element: extracts the trimmed domain, assembles
    ///        integration points via moment fitting, and returns the element.
    /// @param Id      Element id.
    /// @param rBounds Element bounds in global and parametric space.
    /// @return The built element, or std::nullopt if the domain is absent or all IPs are removed.
    [[nodiscard]] std::optional<ElementType> Build(IndexType Id, const ElementBounds& rBounds)
    {
        Timer timer_intersection{};
        auto p_domain = mBrepOperator.get().pGetTrimmedDomain(
            rBounds.global.lower, rBounds.global.upper, mMinVolRatio, mMinNumBoundaryTriangles, mNeglectIfStlIsFlawed
        );
        mElapsedIntersectionTime.fetch_add(timer_intersection.Measure(), std::memory_order_relaxed);

        if (!p_domain) { return std::nullopt; }

        ElementType element(Id, rBounds, std::move(*p_domain));

        Timer timer_fitting{};
        QuadratureTrimmedElement<ElementType>::AssembleIPs(
            element, mPolynomialOrder, mMomentFittingResidual, mEchoLevel
        );
        mElapsedMomentFittingTime.fetch_add(timer_fitting.Measure(), std::memory_order_relaxed);

        if (element.GetIntegrationPoints().empty()) { return std::nullopt; }

        return element;
    }

    /// @brief Returns accumulated time spent in domain extraction.
    /// @details The timing accumulator is updated atomically.
    [[nodiscard]] double ElapsedIntersectionTime() const noexcept
    { return mElapsedIntersectionTime.load(std::memory_order_relaxed); }

    /// @brief Returns accumulated time spent in moment fitting.
    /// @details The timing accumulator is updated atomically.
    [[nodiscard]] double ElapsedMomentFittingTime() const noexcept
    { return mElapsedMomentFittingTime.load(std::memory_order_relaxed); }

    ///@}
private:
    ///@name Private members
    ///@{

    std::reference_wrapper<const BRepOperator> mBrepOperator;
    double mMinVolRatio;
    IndexType mMinNumBoundaryTriangles;
    bool mNeglectIfStlIsFlawed;
    Vector3i mPolynomialOrder;
    double mMomentFittingResidual;
    IndexType mEchoLevel;

    std::atomic<double> mElapsedIntersectionTime{ 0.0 };
    std::atomic<double> mElapsedMomentFittingTime{ 0.0 };
    ///@}
};

/// @class  UntrimmedElementBuilder
/// @brief Encapsulates construction policy for untrimmed elements.
template<
    concepts::IntegrationPoint TIntegrationPointType,
    concepts::BoundaryIntegrationPoint TBoundaryIntegrationPointType>
class UntrimmedElementBuilder
{
public:
    ///@name Type definitions
    ///@{

    using IntegrationPointType = TIntegrationPointType;
    using BoundaryIntegrationPointType = TBoundaryIntegrationPointType;
    using ElementType = UntrimmedElement<IntegrationPointType, BoundaryIntegrationPointType>;
    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;
	using BackgroundGridType = BackgroundGrid<IntegrationPointType, BoundaryIntegrationPointType>; 
	using ElementFilterType = BackgroundGridType::ElementFilter;
	static constexpr ElementFilterType Builds = ElementFilterType::untrimmed;

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructs an UntrimmedElementBuilder from the main settings dictionary.
    /// @param rSettings  The full settings dictionary. Reads non-trimmed quadrature
    ///                   and background-grid sub-dictionaries.
    explicit UntrimmedElementBuilder(const MainDictionaryType& rSettings)
        : mPolynomialOrder(
              rSettings[MainSettings::background_grid_settings].GetRequiredValue<Vector3i>(
                  BackgroundGridSettings::polynomial_order
              )
          ),
          mIntegrationMethod(
              rSettings[MainSettings::non_trimmed_quadrature_rule_settings].GetRequiredValue<IntegrationMethod>(
                  NonTrimmedQuadratureRuleSettings::integration_method
              )
          ),
          mUsesGgqRule(static_cast<int>(mIntegrationMethod) >= 3)
    {}

    ///@}
    ///@name Operations
    ///@{

    /// @brief Builds an untrimmed element and optionally assembles standard Gauss IPs.
    /// @details No-op IP assembly when GGQ is active; GGQ points are assembled later
    ///          outside the element loop via QuadratureMultipleElements.
    /// @param Id      Element id.
    /// @param rBounds Element bounds in global and parametric space.
    /// @return The built element (always valid).
    [[nodiscard]] std::optional<ElementType> Build(IndexType Id, const ElementBounds& rBounds)
    {
        ElementType element(Id, rBounds);
        if (!mUsesGgqRule) {
            QuadratureSingleElement<ElementType>::AssembleIPs(element, mPolynomialOrder, mIntegrationMethod);
        }
        return element;
    }

    /// @brief Returns true if GGQ rules are in use (IPs assembled outside the loop).
    [[nodiscard]] bool UsesGgqRule() const noexcept
    { return mUsesGgqRule; }

    /// @brief Returns the polynomial order stored in this builder.
    [[nodiscard]] const Vector3i& PolynomialOrder() const noexcept
    { return mPolynomialOrder; }

    /// @brief Returns the integration method stored in this builder.
    [[nodiscard]] IntegrationMethod GetIntegrationMethod() const noexcept
    { return mIntegrationMethod; }

    ///@}
private:
    ///@name Private members
    ///@{

    Vector3i mPolynomialOrder;
    IntegrationMethod mIntegrationMethod;
    bool mUsesGgqRule;
    ///@}
};
///@}

}// namespace queso
