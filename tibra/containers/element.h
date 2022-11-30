// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef ELEMENT_INCLUDE_H
#define ELEMENT_INCLUDE_H

// External includes
#include <stdexcept>
#include <memory>

// Project includes
#include "containers/integration_point.h"
#include "utilities/parameters.h"
#include "utilities/mapping_utilities.h"
#include "embedding/trimmed_domain_base.h"

///@name TIBRA Classes
///@{

/**
 * @class  Element
 * @author Manuel Messmer
 * @brief  Element/Knot Spans. Stores quadrature points.
*/
class Element
{

public:

    ///@name Type Defintitions
    ///@{
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef std::array<double,3> PointType;
    typedef std::vector<IntegrationPoint> IntegrationPointVectorType;
    typedef std::vector<std::array<double, 2>> IntegrationPoint1DVectorType;
    typedef std::unique_ptr<TrimmedDomainBase> TrimmedDomainPtrType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@param ElementId UniqueID
    ///@param rLowerBoundParam LowerBound of Element in parametric space.
    ///@param rUpperBoundParam LowerBound of Element in parametric space.
    ///@param rParam TIBRA Parameters.
    Element(IndexType ElementId, const PointType& rLowerBoundParam, const PointType& rUpperBoundParam, const Parameters& rParam) :
        mElementId(ElementId), mLocalLowerPoint(rLowerBoundParam), mLocalUpperPoint(rUpperBoundParam), mParameters(rParam)
    {
        mIsTrimmed = false;
        mTrimmedDomainSetFlag = false;
    }

    // Delete copy constructor
    Element(Element const& rOther) = delete;

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
    /// @return
    const PointType& GetUpperBoundParam() const {
        return mLocalUpperPoint;
    }

    const PointType& GetLowerBoundParam() const {
        return mLocalLowerPoint;
    }

    PointType GetGlobalUpperPoint() const {
        return MappingUtilities::FromLocalToGlobalSpace(mLocalUpperPoint, mParameters.PointA(), mParameters.PointB());
    }

    PointType GetGlobalLowerPoint() const{
        return MappingUtilities::FromLocalToGlobalSpace(mLocalLowerPoint, mParameters.PointA(), mParameters.PointB());
    }

    IntegrationPoint1DVectorType& IntegrationPoints1D(int i){
        if(i ==0)
            return mIntegrationPointsX;
        if(i==1)
            return mIntegrationPointsY;

        return mIntegrationPointsZ;
    }

    IntegrationPoint1DVectorType& IntegrationPointsY(){
        return mIntegrationPointsY;
    }
    IntegrationPoint1DVectorType& IntegrationPointsZ(){
        return mIntegrationPointsZ;
    }

    void pSetTrimmedDomain(TrimmedDomainPtrType& pTrimmedDomain ){
        mpTrimmedDomain = std::move(pTrimmedDomain);
    }

    // Return raw Ptr!!
    const TrimmedDomainBase* const pGetTrimmedDomain() const {
        if( !IsTrimmed() ){
            throw  std::runtime_error("Element :: pGetTrimmedDomain :: Element is not Trimmed." );
        }
        if( !mpTrimmedDomain ){
            throw  std::runtime_error("Element :: pGetTrimmedDomain :: Trimmed Domain Pointer has not been set." );
        }
        return mpTrimmedDomain.get();
    }

    void ClearTrimmedDomain(){
        mpTrimmedDomain = nullptr;
    }

    void SetNeighbourCoefficient(double value, int direction){
        mNumberOfNeighbours[direction] = value;
    }

    double NeighbourCoefficient(){
        return mNumberOfNeighbours[0]*mNumberOfNeighbours[1]*mNumberOfNeighbours[2];
    }

    void SetVisited(bool value){
        mIsVisited = value;
    }

    bool IsVisited(){
        return mIsVisited;
    }

private:
    IntegrationPointVectorType mIntegrationPoints;

    const Parameters& mParameters;
    PointType mLocalUpperPoint;
    PointType mLocalLowerPoint;

    IntegrationPoint1DVectorType mIntegrationPointsX;
    IntegrationPoint1DVectorType mIntegrationPointsY;
    IntegrationPoint1DVectorType mIntegrationPointsZ;

    TrimmedDomainPtrType mpTrimmedDomain = nullptr;
    bool mTrimmedDomainSetFlag = false;

    std::size_t mElementId;
    bool mIsTrimmed;

    std::array<double, 3> mNumberOfNeighbours{};
    bool mIsVisited{};
};

#endif // ELEMENT_INCLUDE_H