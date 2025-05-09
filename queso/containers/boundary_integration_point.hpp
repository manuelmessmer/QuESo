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

#ifndef BOUNDARY_INTEGRATION_POINT_INCLUDE_H
#define BOUNDARY_INTEGRATION_POINT_INCLUDE_H

//// Project includes
#include "queso/containers/integration_point.hpp"

namespace queso {

///@name QuESo Classes
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
    BoundaryIntegrationPoint() : IntegrationPoint(), mNormal{0.0, 0.0, 0.0}
    {}

    /// 3D Constructor
    BoundaryIntegrationPoint(double X, double Y, double Z, double Weigth, const Vector3d& rNormal) :
        IntegrationPoint(X,Y,Z, Weigth), mNormal(rNormal)
    {
    }

    /// Destructor
    ~BoundaryIntegrationPoint() = default;

    /// Copy Constructor
    BoundaryIntegrationPoint(BoundaryIntegrationPoint const& rOther) = default;
    /// Assignement Operator
    BoundaryIntegrationPoint& operator=(BoundaryIntegrationPoint const& rOther) = default;
    /// Move constructor
    BoundaryIntegrationPoint(BoundaryIntegrationPoint&& rOther) noexcept = default;
    /// Move assignement operator
    BoundaryIntegrationPoint& operator=(BoundaryIntegrationPoint&& rOther) noexcept = default;


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
///@} // End QuESo classes

} // End namespace queso

#endif // BOUNDARY_INTEGRATION_POINT_INCLUDE_H