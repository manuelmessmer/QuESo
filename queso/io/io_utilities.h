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

    /// @brief Write points to VTK.
    /// @tparam Type
    /// @param pPoints
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    template<typename Type>
    static void WritePointsToVTK(const std::vector<Type>& rPoints,
                                 const std::string& rFilename,
                                 bool Binary);

private:

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

    ////// Some Helper functions //////

    template<typename T>
    static void SwapEnd(T& var)
    {
        char* varArray = reinterpret_cast<char*>(&var);
        for(long i = 0; i < static_cast<long>(sizeof(var)/2); ++i)
        std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
    }

    template<typename T>
    static void WriteBinary(std::ofstream& stream, T& var){
        SwapEnd(var);
        stream.write(reinterpret_cast<char*>(&var), sizeof(T));
    }

  ///@}
}; // End class IO
///@} End QuESo Classes

} // End namespace queso

#include "queso/io/io_utilities.tpp"

#endif // IO_UTILTIES_H