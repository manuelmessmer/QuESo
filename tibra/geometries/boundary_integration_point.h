// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BOUNDARY_INTEGRATION_POINT_INCLUDE_H
#define BOUNDARY_INTEGRATION_POINT_INCLUDE_H

#include "geometries/integration_point.h"

class BoundaryIntegrationPoint : public IntegrationPoint
{
public:
    // Default constructor
    BoundaryIntegrationPoint() : IntegrationPoint()
    {}

    // 3D Constructor
    BoundaryIntegrationPoint(double x, double y, double z, double weigth_, const std::array<double,3>& rNormal) :
        IntegrationPoint(x,y,z, weigth_), mNormal(rNormal)
    {
    }

    // Assignement operator
    BoundaryIntegrationPoint& operator=(const BoundaryIntegrationPoint& rOther)
    {
        IntegrationPoint::operator=(rOther);
        if( this != &rOther) {
            mNormal[0] = rOther.mNormal[0];
            mNormal[1] = rOther.mNormal[1];
            mNormal[2] = rOther.mNormal[2];
        }
        return *this;
    }

    const std::array<double,3>& Normal(){
        return mNormal;
    }

private:
    std::array<double, 3> mNormal;
};

#endif