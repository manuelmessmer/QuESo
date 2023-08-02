// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BREP_OPERATOR_FACTORY_INCLUDE_H
#define BREP_OPERATOR_FACTORY_INCLUDE_H

#include "embedding/brep_operator.h"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  BRepOperatorFactory
 * @author Manuel Messmer
 * @brief  Provides Functions to Instantiate BRepOperator.
*/
class BRepOperatorFactory {
public:
    ///@name Type Definition
    ///@{

    typedef BRepOperator CurrentBRepOperator;

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns Ptr to new BRepOperator (Unique)
    static Unique<CurrentBRepOperator> New(const TriangleMesh& rTriangleMesh, const Parameters& rParameters){
        return MakeUnique<CurrentBRepOperator>(rTriangleMesh, rParameters);
    }

    ///@}
}; // End BRepOperatorFactory class
///@} // End QuESo classes

} // End namespace queso

#endif