// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef INTEGRATION_POINT_INCLUDE_H
#define INTEGRATION_POINT_INCLUDE_H

#include "containers/point.h"

class IntegrationPoint : public PointType
{
public:
    // Default constructor
    IntegrationPoint() : PointType()
    {}

    // 2D Constructor
    IntegrationPoint(double x, double y, double weigth_) :
        PointType(x,y,0.0), mWeight(weigth_)
    {
        mActiveFlag = true;
    }

    // 3D Constructor
    IntegrationPoint(double x, double y, double z, double weigth_) :
        PointType(x,y,z), mWeight(weigth_)
    {
        mActiveFlag = true;
    }

    // Assignement operator
    IntegrationPoint& operator=(const IntegrationPoint& rOther)
    {
        PointType::operator=(rOther);
        if( this != &rOther) {
            mWeight = rOther.mWeight;
            mActiveFlag = rOther.mActiveFlag;

        }
        return *this;
    }

    double GetWeight() const{
        return mWeight;
    }

    void SetWeight(double weigth_){
        mWeight = weigth_;
    }

    void SetActive(bool value){
        mActiveFlag = value;
    }

    bool IsActive(){
        return mActiveFlag;
    }

private:
    double mWeight;
    bool mActiveFlag;
};

#endif // INTEGRATION_POINT_INCLUDE_H