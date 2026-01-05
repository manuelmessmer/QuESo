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
#include "queso/utilities/math_utilities.hpp"

namespace queso {

// Provides operations two map points between spaces.
namespace Mapping {

    /// @brief Maps point from global to parametric space.
    /// @param rGlobalCoord Point to map.
    /// @param rBoundsXYZ physical bounds of background grid.
    /// @param rBoundsUVW parametric bounds of background grid.
    /// @return PointType
    inline PointType PointFromGlobalToParam( const PointType& rGlobalCoord, const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW){
        const auto delta_xyz = Math::Subtract( rBoundsXYZ.second, rBoundsXYZ.first );
        const auto delta_uvw = Math::Subtract( rBoundsUVW.second, rBoundsUVW.first );

        return PointType{ ( (rGlobalCoord[0] - rBoundsXYZ.first[0]) / std::abs(delta_xyz[0]) * std::abs(delta_uvw[0]) ) + rBoundsUVW.first[0],
                          ( (rGlobalCoord[1] - rBoundsXYZ.first[1]) / std::abs(delta_xyz[1]) * std::abs(delta_uvw[1]) ) + rBoundsUVW.first[1],
                          ( (rGlobalCoord[2] - rBoundsXYZ.first[2]) / std::abs(delta_xyz[2]) * std::abs(delta_uvw[2]) ) + rBoundsUVW.first[2] };
    }

    /// @brief Maps point from parametric to global space.
    /// @param rLocalCoord Point to map.
    /// @param rBoundsXYZ physical bounds of background grid.
    /// @param rBoundsUVW parametric bounds of background grid.
    /// @return PointType.
    inline PointType PointFromParamToGlobal( const PointType& rLocalCoord, const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW){
        const auto delta_xyz = Math::Subtract( rBoundsXYZ.second,  rBoundsXYZ.first );
        const auto delta_uvw = Math::Subtract( rBoundsUVW.second,  rBoundsUVW.first );

        return PointType{ ( (rLocalCoord[0] - rBoundsUVW.first[0]) / std::abs(delta_uvw[0]) * std::abs(delta_xyz[0]) ) + rBoundsXYZ.first[0],
                          ( (rLocalCoord[1] - rBoundsUVW.first[1]) / std::abs(delta_uvw[1]) * std::abs(delta_xyz[1]) ) + rBoundsXYZ.first[1],
                          ( (rLocalCoord[2] - rBoundsUVW.first[2]) / std::abs(delta_uvw[2]) * std::abs(delta_xyz[2]) ) + rBoundsXYZ.first[2] };
    }

    ///@}

} // End namespace Mapping

} // End namespace queso

#endif
