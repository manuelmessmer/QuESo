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

#ifndef MAPPING_UTILITIES_INCLUDE_H
#define MAPPING_UTILITIES_INCLUDE_H

//// Project includes
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

///
/**
 * @class  Mapping
 * @author Manuel Messmer
 * @brief  Provides operations two map between spaces. Static interface.
*/
class Mapping {

public:
    ///@name Public Operations
    ///@{

    /// @brief Maps point from global to parametric space.
    /// @param rGlobalCoord Point to map.
    /// @param rBoundsXYZ physical bounds of background mesh.
    /// @param rBoundsUVW parametric bounds of background mesh.
    /// @return PointType
    static PointType PointFromGlobalToParam( const PointType& rGlobalCoord, const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW);

    /// @brief Maps point from parametric to global space.
    /// @param rLocalCoord Point to map.
    /// @param rBoundsXYZ of background mesh.
    /// @param rBoundsUVW of background mesh.
    /// @return PointType.
    static PointType PointFromParamToGlobal( const PointType& rLocalCoord, const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW);

    ///@}
}; // End class Mapping.
///@} End queso classes.

} // End namespace queso

#endif // End VOXEL_INDEXING_UTILITIES_INCLUDE_H

