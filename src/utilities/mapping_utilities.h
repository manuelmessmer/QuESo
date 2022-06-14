
// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MAPPING_UTILITIES_INCLUDE_H
#define MAPPING_UTILITIES_INCLUDE_H

#include <array>

namespace MappingUtilities
{

std::array<double,3> FromGlobalToLocalSpace(  const std::array<double,3>& rGlobalCoord,
                                              const std::array<double,3>& rLowerPoint,
                                              const std::array<double,3>& rUpperPoint);

std::array<double,3> FromLocalToGlobalSpace( const std::array<double,3>& rLocalCoord,
                                             const std::array<double,3>& rLowerPoint,
                                             const std::array<double,3>& rUpperPoint);

} // end namespace utilities

#endif