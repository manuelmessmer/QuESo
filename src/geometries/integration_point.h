// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef INTEGRATION_POINT_INCLUDE_H
#define INTEGRATION_POINT_INCLUDE_H

#include <memory>
#include "geometries/point3.h"

class IntegrationPoint : public Point3
{
public:
    typedef std::shared_ptr<IntegrationPoint> Pointer;

    // Default constructor
    IntegrationPoint() : Point3()
    {}

    // 2D Constructor
    IntegrationPoint(double x, double y, double weigth_) :
        Point3(x,y,0.0), mWeight(weigth_)
    {
        mActiveFlag = true;
    }

    // 3D Constructor
    IntegrationPoint(double x, double y, double z, double weigth_) :
        Point3(x,y,z), mWeight(weigth_)
    {
        mActiveFlag = true;
    }

    // Assignement operator
    IntegrationPoint& operator=(const IntegrationPoint& rOther)
    {
        Point3::operator=(rOther);
        if( this != &rOther) {
            mWeight = rOther.mWeight;
            mActiveFlag = rOther.mActiveFlag;

        }
        return *this;
    }

    double GetWeight(){
        return mWeight;
    }
    const double GetWeightConst() const {
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