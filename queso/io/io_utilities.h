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

///@name QuESo classes
///@{

/// @class  IO
/// @author Manuel Messmer
/// @brief  Provides methods to read and write data from an to files.
///         Supports STL and VTK file extensions.
class IO {
public:
    ///@}
    ///@name Type definitions
    ///@{

    enum class EncodingType {binary, ascii};

    ///@name Operations
    ///@{

    /*--- Read operations ---*/

    /// @brief Reads TriangleMesh from STL file.
    ///        Checks if file is in binary or ascii format and calls either
    ///        ReadMeshFromSTL_Ascii() or ReadMeshFromSTL_Binary().
    /// @param rTriangleMesh
    /// @param rFilename
    static void ReadMeshFromSTL(TriangleMeshInterface& rTriangleMesh,
                                const std::string& rFilename);

    /*--- Write operations ---*/

    ///@brief Writes dictionary to JSON file.
    ///@tparam TDictType
    ///@param rDictionary
    ///@param rFilename
    template<typename TDictType>
    static void WriteDictionaryToJSON(const TDictType& rDictionary,
                                      const std::string& rFilename);

    /// @brief Writes TriangleMesh to VTK file.
    /// @param rTriangleMesh
    /// @param rFilename
    /// @param Encoding Options: {binary, ascii}.
    static void WriteMeshToVTK(const TriangleMeshInterface& rTriangleMesh,
                               const std::string& rFilename,
                               EncodingType Encoding);

    /// @brief Writes TriangleMesh to STL file.
    /// @param rTriangleMesh
    /// @param rFilename
    /// @param Encoding Options: {binary, ascii}.
    static void WriteMeshToSTL( const TriangleMeshInterface& rTriangleMesh,
                                const std::string& rFilename,
                                EncodingType Encoding);

    /// @brief Writes triangle mesh associated with given condition to STL file.
    /// @tparam TElementType
    /// @param rCondition
    /// @param rFilename
    /// @param Encoding Options: {binary, ascii}.
    template<typename TElementType>
    static void WriteConditionToSTL(const Condition<TElementType>& rCondition,
                                    const std::string& rFilename,
                                    EncodingType Encoding);

    /// @brief Writes elements to VTK file.
    /// @tparam TElementType
    /// @param rBackgroundGrid
    /// @param rFilename
    /// @param Encoding Options: {binary, ascii}.
    template<typename TElementType>
    static void WriteElementsToVTK( const BackgroundGrid<TElementType>& rBackgroundGrid,
                                    const std::string& rFilename,
                                    EncodingType Encoding);

    /// @brief Write points to VTK file.
    /// @tparam TElementType
    /// @param rBackgroundGrid
    /// @param rFilename
    /// @param Encoding Options: {binary, ascii}.
    template<typename TElementType>
    static void WritePointsToVTK(const BackgroundGrid<TElementType>& rBackgroundGrid,
                                 const std::string& rFilename,
                                 EncodingType Encoding);
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

    ///@brief Returns the encopding type of the given file in Ascii-format.
    ///@param rFilename
    ///@return EncodingType: Options {binary, ascii}.
    static EncodingType GetEncodingType(const std::string& rFilename);

    ///@brief  Reads triangle mesh from STL file in Ascii-format.
    ///@param rTriangleMesh
    ///@param rFilename
    ///@see ReadMeshFromSTL_Binary().
    static void ReadMeshFromSTL_Ascii(TriangleMeshInterface& rTriangleMesh,
                                      const std::string& rFilename);

    ///@brief  Reads TriangleMesh from STL file in Binary-format.
    ///@param rTriangleMesh
    ///@param rFilename
    ///@see ReadMeshFromSTL_Ascii().
    static void ReadMeshFromSTL_Binary(TriangleMeshInterface& rTriangleMesh,
                                       const std::string& rFilename);

    /// @brief Helper function to get vertices of hexahedron defined by lower and upper bounds.
    /// @param rLowerBound
    /// @param rUpperBound
    /// @return std::array<PointType, 8>
    static std::array<PointType, 8> GetHexahedronVertices(const PointType& rLowerBound, const PointType& rUpperBound) {
        return {{
            {{rLowerBound[0], rLowerBound[1], rLowerBound[2]}},
            {{rUpperBound[0], rLowerBound[1], rLowerBound[2]}},
            {{rUpperBound[0], rUpperBound[1], rLowerBound[2]}},
            {{rLowerBound[0], rUpperBound[1], rLowerBound[2]}},
            {{rLowerBound[0], rLowerBound[1], rUpperBound[2]}},
            {{rUpperBound[0], rLowerBound[1], rUpperBound[2]}},
            {{rUpperBound[0], rUpperBound[1], rUpperBound[2]}},
            {{rLowerBound[0], rUpperBound[1], rUpperBound[2]}} }};
    }

    ///@}
}; // End class IO
///@} End QuESo Classes

} // End namespace queso

// Include template definitions.
#include "queso/io/io_utilities.tpp"

#endif // IO_UTILTIES_H