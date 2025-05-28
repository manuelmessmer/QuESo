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

#ifndef ELEMENT_INCLUDE_H
#define ELEMENT_INCLUDE_H

//// STL includes

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/containers/data_set.hpp"
#include "queso/embedding/trimmed_domain.h"
#include "queso/utilities/mapping_utilities.hpp"

namespace queso {

/* --- Keys for Element --- */
// Values
QuESo_DEFINE_KEY_SET(ElementValues, MainValuesTypeTag,
    QuESo_KEY_LIST(is_visited, neighbor_coefficient) );
QuESo_DEFINE_KEY_TO_VALUE(ElementValues, is_visited, MainValuesTypeTag, bool);
QuESo_DEFINE_KEY_TO_VALUE(ElementValues, neighbor_coefficient, MainValuesTypeTag, double);
QuESo_REGISTER_KEY_SET(ElementValues, MainValuesTypeTag,
    QuESo_KEY(ElementValues::is_visited),
    QuESo_KEY(ElementValues::neighbor_coefficient)
);

///@name QuESo Classes
///@{

/**
 * @class  Element
 * @author Manuel Messmer
 * @brief  Container to store data on element level.
 *         Is is defined by a simpe bounding box in physical and in parametric space.
 *         Stores quadrature points and trimmed domain (if element is trimmed).
 * @tparam TIntegrationPointType
 * @tparam TBoundaryIntegrationPointType
 * @see    background_grid.hpp
*/
template<typename TIntegrationPointType, typename TBoundaryIntegrationPointType>
class Element
{
public:

    ///@name Type Defintitions
    ///@{

    /// The following point definitions are the general definitions for IntegrationPoint
    /// and BoundaryIntegrationPoint used everywhere else.
    using IntegrationPointType = TIntegrationPointType;
    using BoundaryIntegrationPointType = TBoundaryIntegrationPointType;

    using IntegrationPointVectorType = std::vector<IntegrationPointType>;
    using BoundaryIntegrationPointVectorType = std::vector<BoundaryIntegrationPointType>;
    using TrimmedDomainPtrType = Unique<TrimmedDomain>;

    using DataSetType = DataSet<key::MainValuesTypeTag>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@param ElementId UniqueID
    ///@param rBoundXYZ Bounds of Element in physical space.
    ///@param rBoundUVW Bounds of Element in parametric space.
    Element(IndexType ElementId, const BoundingBoxType& rBoundXYZ, const BoundingBoxType& rBoundUVW) :
        mId(ElementId), mBoundsXYZ(rBoundXYZ), mBoundsUVW(rBoundUVW), mpTrimmedDomain(nullptr),
        mpValues(MakeUnique<DataSetType>(DataSetType::KeySetInfoTypeTag<key::detail::ElementValuesMainValuesTypeTagKeySetInfo>{}))
    {
        mpValues->SetValue(ElementValues::is_visited, false);
        mpValues->SetValue(ElementValues::neighbor_coefficient, 1.0);
    }

    /// Destructor
    ~Element() = default;
    /// Copy constructor
    Element(const Element& rOther) = delete;
    /// Assignement operator
    Element& operator=(const Element& rOther) = delete;
    /// Move constructor
    Element(Element&& rOther) noexcept = default;
    /// Move assignement operator
    Element& operator=(Element&& rOther) noexcept = default;

    ///@}
    ///@name Operations
    ///@{


    /// @brief Return Id of this element.
    /// @return IndexType
    const IndexType GetId() const {
        return mId;
    }

    /// @brief Get bounds of element in physical/global coordinates.
    /// @return BoundingBoxType
    const BoundingBoxType& GetBoundsXYZ() const {
        return mBoundsXYZ;
    }

    /// @brief Get bounds of element in parametric coordinates.
    /// @return BoundingBoxType
    const BoundingBoxType& GetBoundsUVW() const{
        return mBoundsUVW;
    }

    /// @brief Returns Vector of integration points. (non-const)
    /// @return IntegrationPointVectorType&
    IntegrationPointVectorType& GetIntegrationPoints() {
        return mIntegrationPoints;
    }

    /// @brief Returns const& Vector of integration points. (const)
    /// @return const IntegrationPointVectorType&
    const IntegrationPointVectorType& GetIntegrationPoints() const{
        return mIntegrationPoints;
    }


    /// @brief Map point from global space to parametric space.
    /// @param rGlobalCoord
    /// @return PointType
    PointType PointFromGlobalToParam( const PointType& rGlobalCoord ) const {
        return Mapping::PointFromGlobalToParam(rGlobalCoord, mBoundsXYZ, mBoundsUVW);
    }

    /// @brief Map point from parametric space to global space.
    /// @param rLocalCoord
    /// @return PointType
    PointType PointFromParamToGlobal( const PointType& rLocalCoord ) const {
        return Mapping::PointFromParamToGlobal(rLocalCoord, mBoundsXYZ, mBoundsUVW);
    }

    /// @brief Returns determinant of jacobian.
    /// @return double.
    double DetJ() const {
        const auto detla_xyz = Math::Subtract( mBoundsXYZ.second, mBoundsXYZ.first );
        const auto detla_uvw = Math::Subtract( mBoundsUVW.second, mBoundsUVW.first );
        return (detla_xyz[0]*detla_xyz[1]*detla_xyz[2]) / (detla_uvw[0]*detla_uvw[1]*detla_uvw[2]);
    }

    /// @brief Set trimmed domain of element.
    /// @param pTrimmedDomain Ptr (Unique) to new trimmed domain.
    void pSetTrimmedDomain(TrimmedDomainPtrType& pTrimmedDomain ){
        mpTrimmedDomain = std::move(pTrimmedDomain);
    }

    /// @brief Returns true if element is trimmed.
    /// @return bool
    bool IsTrimmed() const {
        return mpTrimmedDomain != nullptr;
    }

    /// @brief Get ptr to trimmed domain of element.
    /// @note Return raw ptr. No transfer of ownership. Element owns trimmed domain.
    /// @return const TrimmedDomain*
    const TrimmedDomain* const pGetTrimmedDomain() const {
        if( !IsTrimmed() ){
            QuESo_ERROR << "Element is not Trimmed.\n";
        }
        return mpTrimmedDomain.get();
    }

    /// @brief Sets the value in the element dataset associated with the given key.
    /// @tparam TKeyType
    /// @tparam TValueType
    /// @param rQueryKey
    /// @param rNewValue
    template<typename TKeyType, typename TValueType>
    void SetValue(const TKeyType& rQueryKey, const TValueType& rNewValue) noexcept(NOTDEBUG) {
        mpValues->SetValue(rQueryKey, rNewValue);
    }

    /// @brief Returns the value in the element dataset associated with the given key.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return const TValueType&
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValue(const TKeyType& rQueryKey) const noexcept(NOTDEBUG) {
        // We can use the fast version. The values are always set, see Constructor.
        return mpValues->GetValueFast<TValueType>(rQueryKey);
    }
    ///@}
private:

    ///@name Private member variables
    ///@{

    const IndexType mId;
    const BoundingBoxType mBoundsXYZ;
    const BoundingBoxType mBoundsUVW;

    IntegrationPointVectorType mIntegrationPoints;
    TrimmedDomainPtrType mpTrimmedDomain;

    Unique<DataSetType> mpValues;
    ///@}
}; // End class Element
///@} // QuESo classes

} // End namespace queso

#endif // ELEMENT_INCLUDE_H