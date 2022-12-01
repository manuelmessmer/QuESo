
// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef PARAMETERS_INCLUDE_H
#define PARAMETERS_INCLUDE_H

// STL includes
#include <stdexcept>
#include <array>

// Project includes
#include "containers/point_types.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace tibra {

class Parameters
{

public:

    // Constructor
    Parameters(PointType PointA, PointType PointB, Vector3i NumberOfElements, Vector3i Order,
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

    const PointType& PointA() const {
        return mPointA;
    }
    const PointType& PointB() const {
        return mPointB;
    }
    const Vector3i& NumberOfElements() const {
        return mNumberOfElements;
    }
    const Vector3i& Order() const {
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
    IntegrationPointFactory1D::IntegrationMethod IntegrationMethod() const {
        // Todo: Make this at Init
        if( mIntegrationMethod == "Gauss" )
            return IntegrationPointFactory1D::IntegrationMethod::Gauss;
        else if( mIntegrationMethod == "ReducedGauss1" )
            return IntegrationPointFactory1D::IntegrationMethod::ReducedGauss1;
        else if( mIntegrationMethod == "ReducedGauss2" )
            return IntegrationPointFactory1D::IntegrationMethod::ReducedGauss2;
        else if( mIntegrationMethod == "ReducedExact")
            return IntegrationPointFactory1D::IntegrationMethod::ReducedExact;
        else if( mIntegrationMethod == "ReducedOrder1")
            return IntegrationPointFactory1D::IntegrationMethod::ReducedOrder1;
        else if( mIntegrationMethod == "ReducedOrder2")
            return IntegrationPointFactory1D::IntegrationMethod::ReducedOrder2;
        else
            throw  std::invalid_argument("Parameters: Integration Method: " + mIntegrationMethod + " not available! \n");
    }
    int EchoLevel() const {
        return mEchoLevel;
    }
    void SetUseCustomizedTrimmedPointsPositionFlag(bool value){
        mUseCustomizedTrimmedPointsPositionFlag = value;
    }
    bool UseCustomizedTrimmedPositions() const{
        return mUseCustomizedTrimmedPointsPositionFlag;
    }
    const double GetPointDistributionFactor() const {
        return mPointDistributionFactor;
    }

private:
    const PointType mPointA;
    const PointType mPointB;
    const Vector3i mNumberOfElements;
    const Vector3i mOrder;
    const double mInitialTriangleEdgeLength;
    const int mMinimumNumberOfTriangles;
    const double mMomentFittingResidual;
    const std::string mIntegrationMethod;
    const int mEchoLevel;
    const double mPointDistributionFactor;
    bool mUseCustomizedTrimmedPointsPositionFlag;
};

} // End namespace tibra

#endif

