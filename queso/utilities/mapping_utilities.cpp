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

// STL includes

// Project includes
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/mapping_utilities.h"

namespace queso {

PointType Mapping::PointFromGlobalToParam( const PointType& rGlobalCoord, const BoundingBoxType& rBoundXYZ, const BoundingBoxType& rBoundUVW){
    const auto delta_xyz = Math::Subtract( rBoundXYZ.second, rBoundXYZ.first );
    const auto delta_uvw = Math::Subtract( rBoundUVW.second, rBoundUVW.first );

    return PointType{ ( (rGlobalCoord[0] - rBoundXYZ.first[0]) / std::abs(delta_xyz[0]) * std::abs(delta_uvw[0]) ) + rBoundUVW.first[0],
                      ( (rGlobalCoord[1] - rBoundXYZ.first[1]) / std::abs(delta_xyz[1]) * std::abs(delta_uvw[1]) ) + rBoundUVW.first[1],
                      ( (rGlobalCoord[2] - rBoundXYZ.first[2]) / std::abs(delta_xyz[2]) * std::abs(delta_uvw[2]) ) + rBoundUVW.first[2] };
}

PointType Mapping::PointFromParamToGlobal( const PointType& rLocalCoord, const BoundingBoxType& rBoundXYZ, const BoundingBoxType& rBoundUVW ) {
    const auto delta_xyz = Math::Subtract( rBoundXYZ.second,  rBoundXYZ.first );
    const auto delta_uvw = Math::Subtract( rBoundUVW.second,  rBoundUVW.first );

    return PointType{ ( (rLocalCoord[0] - rBoundUVW.first[0]) / std::abs(delta_uvw[0]) * std::abs(delta_xyz[0]) ) + rBoundXYZ.first[0],
                      ( (rLocalCoord[1] - rBoundUVW.first[1]) / std::abs(delta_uvw[1]) * std::abs(delta_xyz[1]) ) + rBoundXYZ.first[1],
                      ( (rLocalCoord[2] - rBoundUVW.first[2]) / std::abs(delta_uvw[2]) * std::abs(delta_xyz[2]) ) + rBoundXYZ.first[2] };
}

} // End namespace queso