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

#ifndef IO_UTILTIES_H
#define IO_UTILTIES_H

//// STL includes
#include <fstream>
#include <numeric>

//// Project includes
#include "queso/containers/background_grid.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/containers/boundary_integration_point.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  IO
 * @author Manuel Messmer
 * @brief  Provides methods to parse data. Supports STL and VTK files.
**/
class IO {

public:
    ///@name Operations
    ///@{

    /// @brief Write TriangleMeshInterface to VTK-File
    /// @param rTriangleMesh
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    static void WriteMeshToVTK(const TriangleMeshInterface& rTriangleMesh,
                               const std::string& rFilename,
                               bool Binary);

    /// @brief Write TriangleMeshInterface to STL-File.
    /// @param rTriangleMesh
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    static void WriteMeshToSTL( const TriangleMeshInterface& rTriangleMesh,
                                const std::string& rFilename,
                                bool Binary);

    /// @brief Read TriangleMeshInterface from STL.
    /// @param rTriangleMesh
    /// @param rFilename
    static void ReadMeshFromSTL(TriangleMeshInterface& rTriangleMesh,
                                const std::string& rFilename);

    ///@brief Write dictionary to JSON file.
    ///@tparam TDictType
    ///@param rDictionary
    ///@param rFilename
    template<typename TDictType>
    static void WriteDictionaryToJSON(const TDictType& rDictionary, const std::string& rFilename);

    /// @brief Write triangle mesh associated to given condition to STL-file.
    /// @tparam TElementType
    /// @param rCondition
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    template<typename TElementType>
    static void WriteConditionToSTL(const Condition<TElementType>& rCondition,
                                    const std::string& rFilename,
                                    bool Binary);

    /// @brief Write element container to VTK-file.
    /// @tparam TElementType
    /// @param rBackgroundGrid
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    template<typename TElementType>
    static void WriteElementsToVTK( const BackgroundGrid<TElementType>& rBackgroundGrid,
                                    const std::string& rFilename,
                                    bool Binary);

    /// @brief Write points to VTK. Interface for BackgroundGrid.
    /// @tparam TElementType
    /// @param rBackgroundGrid
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    /// @todo Needs to be refactored.
    template<typename TElementType>
    static void WritePointsToVTK(const BackgroundGrid<TElementType>& rBackgroundGrid,
                                 const std::string& rFilename,
                                 bool Binary);

private:
    ///@}
    ///@name Type definitions
    ///@{

    struct PointComparison {
        bool operator() (const PointType& lhs, const PointType& rhs) const {
            const double dx = lhs[0] - rhs[0];
            const double dy = lhs[1] - rhs[1];
            const double dz = lhs[2] - rhs[2];

            if (std::abs(dx) < SNAPTOL) {
                if (std::abs(dy) < SNAPTOL) {
                    // x and y close enough, compare z
                    return dz < -SNAPTOL;
                }
                // x close enough, compare y
                return dy < -SNAPTOL;
            }
            // compare x
            return dx < -SNAPTOL;
        }
    };

    ///@}
    ///@name Private Operations
    ///@{

    ///@brief  Reads triangle mesh from STL in Ascii-format.
    ///@param rTriangleMesh
    ///@param rFilename
    ///@see ReadMeshFromSTL_Binary()
    static void ReadMeshFromSTL_Ascii(TriangleMeshInterface& rTriangleMesh,
                                      const std::string& rFilename);

                                      ///@brief  Reads triangle mesh from STL in Binary-format.
    ///@param rTriangleMesh
    ///@param rFilename
    ///@see ReadMeshFromSTL_Ascii()
    static void ReadMeshFromSTL_Binary(TriangleMeshInterface& rTriangleMesh,
        const std::string& rFilename);

    ///@brief Returns true if given file in in ASCII-format.
    ///@param rFilename
    ///@return bool
    static bool STLIsInASCIIFormat(const std::string& rFilename);

    /// @brief Helper function to get hexahedron vertices from bounds.
    /// @param rMin
    /// @param rMax
    /// @return std::array<PointType, 8>
    static std::array<PointType, 8> GetHexahedronVertices(const PointType& rMin, const PointType& rMax) {
        return {{
            {{rMin[0], rMin[1], rMin[2]}},
            {{rMax[0], rMin[1], rMin[2]}},
            {{rMax[0], rMax[1], rMin[2]}},
            {{rMin[0], rMax[1], rMin[2]}},
            {{rMin[0], rMin[1], rMax[2]}},
            {{rMax[0], rMin[1], rMax[2]}},
            {{rMax[0], rMax[1], rMax[2]}},
            {{rMin[0], rMax[1], rMax[2]}} }};
    }

    ///@}
}; // End class IO
///@} End QuESo Classes

} // End namespace queso

#include "queso/io/io_utilities.tpp"

#endif // IO_UTILTIES_H