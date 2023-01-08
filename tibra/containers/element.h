// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef ELEMENT_INCLUDE_H
#define ELEMENT_INCLUDE_H

//// STL includes
#include <stdexcept>
#include <memory>
//// Project includes
#include "containers/integration_point.h"
#include "embedding/trimmed_domain_base.h"
#include "utilities/parameters.h"
#include "utilities/mapping_utilities.h"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  Element
 * @author Manuel Messmer
 * @brief  Element/Knot Spans. Stores quadrature points and trimmed domain (if element is trimmed).
*/
class Element
{

public:

    ///@name Type Defintitions
    ///@{
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef std::vector<IntegrationPoint> IntegrationPointVectorType;
    typedef std::vector<std::array<double, 2>> IntegrationPoint1DVectorType;
    typedef Unique<TrimmedDomainBase> TrimmedDomainPtrType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@param ElementId UniqueID
    ///@param rLowerBoundParam LowerBound of Element in parametric space.
    ///@param rUpperBoundParam LowerBound of Element in parametric space.
    ///@param rParam TIBRA Parameters.
    Element(IndexType ElementId, const PointType& rLowerBoundParam, const PointType& rUpperBoundParam, const Parameters& rParam) :
        mElementId(ElementId), mLowerBoundParam(rLowerBoundParam), mUpperBoundParam(rUpperBoundParam), mParameters(rParam)
    {
        mIsTrimmed = false;
    }

    // Delete copy constructor
    Element(Element const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Set Element as trimmed.
    /// @param Value
    void SetIsTrimmed(bool Value){
        mIsTrimmed = Value;
    }

    /// @brief Set Id
    /// @param Value
    void SetId(IndexType Value){
        mElementId = Value;
    }

    /// @brief Return Id of this element.
    /// @return IndexType
    const IndexType GetId() const{
        return mElementId;
    }

    /// @brief Returns true if element is trimmed.
    /// @return bool
    bool IsTrimmed() const{
        return mIsTrimmed;
    }

    /// @brief Returns TIBRA parameters
    /// @return const Parameters&
    const Parameters& GetParameters() const {
        return mParameters;
    }

    /// @brief Returns Vector of integration points.
    /// @return IntegrationPointVectorType&
    IntegrationPointVectorType& GetIntegrationPoints() {
        return mIntegrationPoints;
    }

    /// @brief Returns const& Vector of integration points.
    /// @return const IntegrationPointVectorType&
    const IntegrationPointVectorType& GetIntegrationPoints() const{
        return mIntegrationPoints;
    }

    /// @brief Get UpperBound of element in parametric coordinates.
    /// @return const PointType&
    const PointType& GetUpperBoundParam() const {
        return mUpperBoundParam;
    }

    /// @brief Get LowerBound of element in parametric coordinates.
    /// @return const PointType&
    const PointType& GetLowerBoundParam() const {
        return mLowerBoundParam;
    }

    /// @brief Get UpperBound of element in physical/global coordinates.
    /// @return PointType
    PointType GetUpperBound() const {
        return Mapping::ParamToGlobal(mUpperBoundParam, mParameters.LowerBound(), mParameters.UpperBound());
    }

    /// @brief Get LowerBound of element in physical/global coordinates.
    /// @return PointType
    PointType GetLowerBound() const{
        return Mapping::ParamToGlobal(mLowerBoundParam, mParameters.LowerBound(), mParameters.UpperBound());
    }

    /// @brief Returns 1D integration points. Required for assembly of GGQ rules.
    /// @param Dir Space Direction: 0-x, 1-y, 2-z.
    /// @return IntegrationPoint1DVectorType&
    IntegrationPoint1DVectorType& IntegrationPoints1D(IndexType Dir){
        if(Dir==0)
            return mIntegrationPointsX;
        if(Dir==1)
            return mIntegrationPointsY;

        return mIntegrationPointsZ;
    }

    /// @brief Set trimmed domain of element.
    /// @param pTrimmedDomain Ptr (Unique) to new trimmed domain.
    void pSetTrimmedDomain(TrimmedDomainPtrType& pTrimmedDomain ){
        mpTrimmedDomain = std::move(pTrimmedDomain);
    }

    /// @brief Get ptr to trimmed domain of element.
    /// @note Return raw ptr. No transfer of ownership. Element owns trimmed domain.
    /// @return const TrimmedDomainBase*
    const TrimmedDomainBase* const pGetTrimmedDomain() const {
        if( !IsTrimmed() ){
            throw  std::runtime_error("Element :: pGetTrimmedDomain :: Element is not Trimmed." );
        }
        if( !mpTrimmedDomain ){
            throw  std::runtime_error("Element :: pGetTrimmedDomain :: Trimmed Domain Pointer has not been set." );
        }
        return mpTrimmedDomain.get();
    }

    /// @brief Clear trimmed domain of element.
    void ClearTrimmedDomain(){
        mpTrimmedDomain = nullptr;
    }

    /// @brief Set neighbour coefficient. Required for assembly of GGQ rule. See: multiple_elements.h.
    /// @param Value New Value.
    /// @param Direction Space Direction: 0-x, 1-y, 2-z.
    void SetNeighbourCoefficient(double Value, IndexType Direction){
        mNumberOfNeighbours[Direction] = Value;
    }

    /// @brief Get neighbour coeefficient of this element. Required for assembly of GGQ rule. See: multiple_elements.h.
    /// @return double
    double NeighbourCoefficient(){
        return mNumberOfNeighbours[0]*mNumberOfNeighbours[1]*mNumberOfNeighbours[2];
    }

    /// @brief Set Flag.
    /// @param Value
    void SetVisited(bool Value){
        mIsVisited = Value;
    }

    /// @brief Returns Flag. (see. SetVisited(value)).
    /// @return bool
    bool IsVisited(){
        return mIsVisited;
    }

    ///@}
private:

    ///@name Private member variables
    ///@{
    IntegrationPointVectorType mIntegrationPoints{};

    IntegrationPoint1DVectorType mIntegrationPointsX{};
    IntegrationPoint1DVectorType mIntegrationPointsY{};
    IntegrationPoint1DVectorType mIntegrationPointsZ{};

    const Parameters& mParameters;
    PointType mUpperBoundParam{};
    PointType mLowerBoundParam{};

    TrimmedDomainPtrType mpTrimmedDomain = nullptr;

    IndexType mElementId{};
    bool mIsTrimmed{};
    bool mIsVisited{};

    PointType mNumberOfNeighbours{};
    ///@}
}; // End class Element
///@} // TIBRA classes

} // End namespace tibra

#endif // ELEMENT_INCLUDE_H