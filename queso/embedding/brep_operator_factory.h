// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BREP_OPERATOR_FACTORY_INCLUDE_H
#define BREP_OPERATOR_FACTORY_INCLUDE_H

#if defined USE_CGAL
#include "cgal_wrapper/cgal_brep_operator.h"
#else
#include "embedding/brep_operator.h"
#endif

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  BRepOperatorFactory
 * @author Manuel Messmer
 * @brief  Provides Functions to Instantiate BRepOperator
 * @details Used to choose between BRepOperator and cgal::CGALBRepOperator.
*/
class BRepOperatorFactory {
public:
    ///@name Type Definition
    ///@{

    #if defined USE_CGAL
    typedef cgal::CGALBRepOperator CurrentBRepOperator;
    #else
    typedef BRepOperator CurrentBRepOperator;
    #endif

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