
// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef PARAMETERS_INCLUDE_H
#define PARAMETERS_INCLUDE_H

#include <array>
#include "utilities/integration_points/integration_points_factory.h"

typedef std::array<double, 3> DoublePointType;
typedef std::array<int, 3> IntPointType;

class Parameters
{
public:
    Parameters(DoublePointType PointA, DoublePointType PointB, IntPointType NumberOfElements, IntPointType Order,
               double InitialTriangleEdgeLength, int MinimumNumberOfTriangles, double MomentFittingResidual,
               double PointDistributionFactor, std::string IntegrationMethod, int EchoLevel) :
            mPointA(PointA),
            mPointB(PointB),
            mNumberOfElements(NumberOfElements),
            mOrder(Order),
            mInitialTriangleEdgeLength(InitialTriangleEdgeLength),
            mMinimumNumberOfTriangles(MinimumNumberOfTriangles),
            mMomentFittingResidual(MomentFittingResidual),
            mPointDistributionFactor(PointDistributionFactor),
            mIntegrationMethod(IntegrationMethod),
            mEchoLevel(EchoLevel)
    {
        mUseCustomizedTrimmedPointsPositionFlag = false;
    }

    const DoublePointType& PointA() const {
        return mPointA;
    }
    const DoublePointType& PointB() const {
        return mPointB;
    }
    const IntPointType& NumberOfElements() const {
        return mNumberOfElements;
    }
    const IntPointType& Order() const {
        return mOrder;
    }
    double InitialTriangleEdgeLength() const {
        return mInitialTriangleEdgeLength;
    }
    int MinimumNumberOfTriangles() const {
        return mMinimumNumberOfTriangles;
    }
    double MomentFittingResidual() const {
        return mMomentFittingResidual;
    }
    IntegrationPointFactory::IntegrationMethod IntegrationMethod() const {
        if( mIntegrationMethod == "Gauss" )
            return IntegrationPointFactory::IntegrationMethod::Gauss;
        else if( mIntegrationMethod == "ReducedGauss1" )
            return IntegrationPointFactory::IntegrationMethod::ReducedGauss1;
        else if( mIntegrationMethod == "ReducedGauss2" )
            return IntegrationPointFactory::IntegrationMethod::ReducedGauss2;
        else if( mIntegrationMethod == "ReducedExact")
            return IntegrationPointFactory::IntegrationMethod::ReducedExact;
        else if( mIntegrationMethod == "ReducedOrder1")
            return IntegrationPointFactory::IntegrationMethod::ReducedOrder1;
        else if( mIntegrationMethod == "ReducedOrder2")
            return IntegrationPointFactory::IntegrationMethod::ReducedOrder2;
        else
            throw  std::runtime_error("Parameters: Integration Method: " + mIntegrationMethod + " not available! \n");
    }
    int EchoLevel() const {
        return mEchoLevel;
    }
    void SetUseCustomizedTrimmedPointsPositionFlag(bool value){
        mUseCustomizedTrimmedPointsPositionFlag = true;
    }
    bool UseCustomizedTrimmedPositions() const{
        return mUseCustomizedTrimmedPointsPositionFlag;
    }
    const double GetPointDistributionFactor() const {
        return mPointDistributionFactor;
    }

private:
    const DoublePointType mPointA;
    const DoublePointType mPointB;
    const IntPointType mNumberOfElements;
    const IntPointType mOrder;
    const double mInitialTriangleEdgeLength;
    const int mMinimumNumberOfTriangles;
    const double mMomentFittingResidual;
    const std::string mIntegrationMethod;
    const int mEchoLevel;
    const double mPointDistributionFactor;
    bool mUseCustomizedTrimmedPointsPositionFlag;
};

#endif

