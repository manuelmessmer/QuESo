// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef ELEMENT_H
#define ELEMENT_H

// External includes
#include <stdexcept>
#include "omp.h"
#include <memory>

// Project includes
#include "geometries/integration_point.h"
#include "utilities/parameters.h"
#include "utilities/mapping_utilities.h"
#include "embedding/trimmed_domain_base.h"

class Element
{

public:

    // Typedefs
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef std::array<double,3> PointType;
    typedef std::vector<IntegrationPoint> IntegrationPointVectorType;
    typedef std::vector<std::array<double, 2>> IntegrationPoint1DVectorType;
    typedef std::unique_ptr<TrimmedDomainBase> TrimmedDomainPtrType;

    // Constructor
    Element(std::size_t ID, PointType PointLocalA, PointType PointLocalB, const Parameters& rParam) :
        mElementId(ID), mLocalLowerPoint(PointLocalA), mLocalUpperPoint(PointLocalB), mParameters(rParam)
    {
        mIsTrimmed = false;
        mTrimmedDomainSetFlag = false;
    }

    // Delete copy constructor
    Element(Element const& rOther) = delete;

    void SetIsTrimmed(bool value){
        mIsTrimmed = value;
    }

    void SetId(std::size_t value){
        mElementId = value;
    }

    const std::size_t GetId() const{
        return mElementId;
    }

    bool IsTrimmed() const{
        return mIsTrimmed;
    }

    const Parameters& GetParameters() const {
        return mParameters;
    }

    IntegrationPointVectorType& GetIntegrationPointsTrimmed(){
        return mIntegrationPointsTrimmed;
    }
    IntegrationPointVectorType& GetIntegrationPointsInside(){
        return mIntegrationPointsInside;
    }
    IntegrationPointVectorType& GetIntegrationPointsFictitious(){
        return mIntegrationPointsFictitious;
    }

    PointType GetLocalUpperPoint() const {
        return mLocalUpperPoint;
    }

    PointType GetLocalLowerPoint() const {
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
    IntegrationPointVectorType mIntegrationPointsTrimmed;
    IntegrationPointVectorType mIntegrationPointsInside;
    IntegrationPointVectorType mIntegrationPointsFictitious;

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

#endif