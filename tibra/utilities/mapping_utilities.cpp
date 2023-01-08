// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// STL includes
#include <cstdlib>

// Project includes
#include "utilities/mapping_utilities.h"

namespace tibra {

PointType Mapping::GlobalToParam( const PointType& rGlobalCoord,
                                  const PointType& rLowerPoint,
                                  const PointType& rUpperPoint){

    return PointType( (rGlobalCoord[0] - rLowerPoint[0]) / std::abs(rUpperPoint[0] - rLowerPoint[0]),
                      (rGlobalCoord[1] - rLowerPoint[1]) / std::abs(rUpperPoint[1] - rLowerPoint[1]),
                      (rGlobalCoord[2] - rLowerPoint[2]) / std::abs(rUpperPoint[2] - rLowerPoint[2]) );
}

PointType Mapping::ParamToGlobal( const PointType& rLocalCoord,
                                  const PointType& rLowerPoint,
                                  const PointType& rUpperPoint ) {

    return PointType( (rLocalCoord[0] * std::abs(rUpperPoint[0] - rLowerPoint[0]) ) + rLowerPoint[0],
                      (rLocalCoord[1] * std::abs(rUpperPoint[1] - rLowerPoint[1]) ) + rLowerPoint[1],
                      (rLocalCoord[2] * std::abs(rUpperPoint[2] - rLowerPoint[2]) ) + rLowerPoint[2] );
}

} // End namespace tibra