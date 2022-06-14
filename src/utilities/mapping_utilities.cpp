
// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#include "utilities/mapping_utilities.h"

std::array<double,3> MappingUtilities::FromGlobalToLocalSpace( const std::array<double,3>& rGlobalCoord,
                                              const std::array<double,3>& rLowerPoint,
                                              const std::array<double,3>& rUpperPoint){
    std::array<double,3> local_coord;
    local_coord[0] = (rGlobalCoord[0] - rLowerPoint[0]) / std::abs(rUpperPoint[0] - rLowerPoint[0]);
    local_coord[1] = (rGlobalCoord[1] - rLowerPoint[1]) / std::abs(rUpperPoint[1] - rLowerPoint[1]);
    local_coord[2] = (rGlobalCoord[2] - rLowerPoint[2]) / std::abs(rUpperPoint[2] - rLowerPoint[2]);

    return local_coord;
}

std::array<double,3> MappingUtilities::FromLocalToGlobalSpace( const std::array<double,3>& rLocalCoord,
                                              const std::array<double,3>& rLowerPoint,
                                              const std::array<double,3>& rUpperPoint){
    std::array<double,3> global_cood;
    global_cood[0] = (rLocalCoord[0] * std::abs(rUpperPoint[0] - rLowerPoint[0]) ) + rLowerPoint[0];
    global_cood[1] = (rLocalCoord[1] * std::abs(rUpperPoint[1] - rLowerPoint[1]) ) + rLowerPoint[1];
    global_cood[2] = (rLocalCoord[2] * std::abs(rUpperPoint[2] - rLowerPoint[2]) ) + rLowerPoint[2];

    return global_cood;
}