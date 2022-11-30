// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// External includes
#include <cstdlib>

// Project includes
#include "utilities/mapping_utilities.h"

PointType MappingUtilities::FromGlobalToLocalSpace( const PointType& rGlobalCoord,
                                              const PointType& rLowerPoint,
                                              const PointType& rUpperPoint){
    PointType local_coord;
    local_coord[0] = (rGlobalCoord[0] - rLowerPoint[0]) / std::abs(rUpperPoint[0] - rLowerPoint[0]);
    local_coord[1] = (rGlobalCoord[1] - rLowerPoint[1]) / std::abs(rUpperPoint[1] - rLowerPoint[1]);
    local_coord[2] = (rGlobalCoord[2] - rLowerPoint[2]) / std::abs(rUpperPoint[2] - rLowerPoint[2]);

    return local_coord;
}

PointType MappingUtilities::FromLocalToGlobalSpace( const PointType& rLocalCoord,
                                              const PointType& rLowerPoint,
                                              const PointType& rUpperPoint){
    PointType global_cood;
    global_cood[0] = (rLocalCoord[0] * std::abs(rUpperPoint[0] - rLowerPoint[0]) ) + rLowerPoint[0];
    global_cood[1] = (rLocalCoord[1] * std::abs(rUpperPoint[1] - rLowerPoint[1]) ) + rLowerPoint[1];
    global_cood[2] = (rLocalCoord[2] * std::abs(rUpperPoint[2] - rLowerPoint[2]) ) + rLowerPoint[2];

    return global_cood;
}