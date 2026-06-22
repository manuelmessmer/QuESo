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
#include "queso/containers/cell_domain.hpp"
#include "queso/containers/data_set.hpp"
#include "queso/containers/element_base.hpp"
#include "queso/includes/define.hpp"

namespace queso {

QuESo_DEFINE_KEY_SET(ElementValues, MainValuesTypeTag, QuESo_KEY_LIST(is_visited, neighbor_coefficient));
QuESo_DEFINE_KEY_TO_VALUE(ElementValues, is_visited, MainValuesTypeTag, bool, KeyRequirement::required);
QuESo_DEFINE_KEY_TO_VALUE(ElementValues, neighbor_coefficient, MainValuesTypeTag, double, KeyRequirement::required);
QuESo_REGISTER_KEY_SET(
    ElementValues,
    MainValuesTypeTag,
    QuESo_KEY(ElementValues::is_visited),
    QuESo_KEY(ElementValues::neighbor_coefficient)
);

///@name QuESo Classes
///@{

/// @class UntrimmedElement
/// @author Manuel Messmer
/// @brief Element representation for untrimmed background-grid cells.
/// @details Extends ElementBase with a CellDomain and stores auxiliary
/// per-element data used during grid traversal and neighborhood operations.
/// The active domain coincides with the full element bounds.
/// @tparam TIntegrationPoint Volume integration-point type.
/// @tparam TBoundaryIntegrationPoint Boundary integration-point type.
template<typename TIntegrationPoint, typename TBoundaryIntegrationPoint>
class UntrimmedElement : public ElementBase<TIntegrationPoint, TBoundaryIntegrationPoint, CellDomain, false>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ElementBase<TIntegrationPoint, TBoundaryIntegrationPoint, CellDomain, false>;

    using CoreType = typename BaseType::CoreType;
    using ElementViewType = typename BaseType::ElementViewType;
    using IntegrationPointType = typename BaseType::IntegrationPointType;
    using BoundaryIntegrationPointType = typename BaseType::BoundaryIntegrationPointType;
    using IntegrationPointVectorType = typename BaseType::IntegrationPointVectorType;
    using BoundaryIntegrationPointVectorType = typename BaseType::BoundaryIntegrationPointVectorType;
    using DataSetType = DataSet<key::MainValuesTypeTag>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructs an untrimmed element from its components.
    /// @param Id Unique element identifier.
    /// @param rBounds Element bounds in global and parametric space.
    /// @param IntegrationPoints Volume integration points in parametric space.
    UntrimmedElement(IndexType Id, const ElementBounds& rBounds, IntegrationPointVectorType IntegrationPoints = {})
        : UntrimmedElement(CoreType{ Id, rBounds, std::move(IntegrationPoints) })
    {}

    /// @brief Constructs an untrimmed element from a preassembled core.
    /// @param rCore Element core data.
    explicit UntrimmedElement(CoreType&& rCore)
        : BaseType(std::move(rCore), CellDomain{}),
          mDataSet(DataSetType::KeySetInfoTypeTag<key::detail::ElementValuesMainValuesTypeTagKeySetInfo>{})
    {
        mDataSet.SetValue(ElementValues::is_visited, false);
        mDataSet.SetValue(ElementValues::neighbor_coefficient, 1.0);
        mDataSet.CheckRequired();
    }

    ///@}
    ///@name Data Access
    ///@{

    /// @brief Sets a value in the element data set.
    /// @tparam TKeyType Key type.
    /// @tparam TValueType Value type.
    /// @param rQueryKey Data-set key.
    /// @param rNewValue New value.
    template<typename TKeyType, typename TValueType>
    void SetValue(const TKeyType& rQueryKey, const TValueType& rNewValue) noexcept(NOTDEBUG)
    { mDataSet.SetValue(rQueryKey, rNewValue); }

    /// @brief Returns a value stored in the element data set.
    /// @tparam TValueType Expected value type.
    /// @tparam TKeyType Key type.
    /// @param rQueryKey Data-set key.
    /// @return const TValueType&
    template<typename TValueType, typename TKeyType>
    [[nodiscard]] const TValueType& GetValue(const TKeyType& rQueryKey) const noexcept(NOTDEBUG)
    { return mDataSet.GetRequiredValue<TValueType>(rQueryKey); }

    ///@}

private:
    ///@name Private Members
    ///@{

    /// @brief Auxiliary per-element data storage.
    DataSetType mDataSet;

    ///@}
};

///@}

}// namespace queso
