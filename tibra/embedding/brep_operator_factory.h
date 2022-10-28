// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BREP_OPERATOR_FACTORY_INCLUDE_H
#define BREP_OPERATOR_FACTORY_INCLUDE_H

#if defined USE_CGAL
#include "cgal_wrapper/cgal_brep_operator.h"
#endif

///@name TIBRA Classes
///@{

/**
 * @class  BRepOperatorFactory
 * @author Manuel Messmer
 * @brief  Provides Functions to Insntiate BRepOperator
 * @details Used to choose between BRepOperator and cgal::BRepOperator.
*/
class BRepOperatorFactory {
public:
    ///@name Type Definition
    ///@{

    #if defined USE_CGAL
    typedef cgal::BRepOperator CurrentBRepOperator;
    #endif

    ///@}
    ///@name Type Definition
    ///@{

    ///@brief Returns Ptr to new BRepOperator
    static std::unique_ptr<CurrentBRepOperator> New(const TriangleMesh& rTriangleMesh){
        return std::make_unique<CurrentBRepOperator>(rTriangleMesh);
    }

    ///@}
};
///@}

#endif