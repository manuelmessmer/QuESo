// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BREP_OPERATOR_FACTORY_INCLUDE_H
#define BREP_OPERATOR_FACTORY_INCLUDE_H

#if defined USE_CGAL
#include "cgal_wrapper/cgal_brep_operator.h"
#endif

namespace tibra {

///@name TIBRA Classes
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
    #endif

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns Ptr to new BRepOperator (std::unique_ptr)
    static std::unique_ptr<CurrentBRepOperator> New(const TriangleMesh& rTriangleMesh){
        return std::make_unique<CurrentBRepOperator>(rTriangleMesh);
    }

    ///@}
}; // End BRepOperatorFactory class
///@} // End TIBRA classes

} // End namespace tibra

#endif