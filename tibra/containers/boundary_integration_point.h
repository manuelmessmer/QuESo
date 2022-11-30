// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BOUNDARY_INTEGRATION_POINT_INCLUDE_H
#define BOUNDARY_INTEGRATION_POINT_INCLUDE_H

#include "containers/integration_point.h"

///@name TIBRA Classes
///@{

/**
 * @class  BoundaryIntegrationPoint
 * @author Manuel Messmer
 * @brief Derives from IntegrationPoint. Contains also Normal Vector.
*/
class BoundaryIntegrationPoint : public IntegrationPoint
{
public:

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    BoundaryIntegrationPoint() : IntegrationPoint()
    {}

    /// 3D Constructor
    /// @param x
    /// @param y
    /// @param z
    /// @param Weigth
    /// @param rNormal
    BoundaryIntegrationPoint(double X, double Y, double Z, double Weigth, const Vector3d& rNormal) :
        IntegrationPoint(X,Y,Z, Weigth), mNormal(rNormal)
    {
    }

    /// Assignement operator
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

    ///@}
    ///@name Operations
    ///@{

    /// @brief Get Normal vector
    /// @return
    const Vector3d& Normal() const{
        return mNormal;
    }

    ///@}

private:
    ///@name Private Member Variables
    ///@{
    Vector3d mNormal;
    ///@}
}; // End class BoundaryIntegrationPoint
///@}

#endif // BOUNDARY_INTEGRATION_POINT_INCLUDE_H