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
class IO{

public:
    ///@name Operations
    ///@{

    /// @brief Write TriangleMeshInterface to VTK-File
    /// @param rTriangleMesh
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    static void WriteMeshToVTK(const TriangleMeshInterface& rTriangleMesh,
                                const std::string& rFilename,
                                const bool Binary);

    /// @brief Write TriangleMeshInterface to STL-File.
    /// @param rTriangleMesh
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    static void WriteMeshToSTL( const TriangleMeshInterface& rTriangleMesh,
                                const std::string& rFilename,
                                const bool Binary);

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
    static void WriteDictionaryToJSON(const TDictType& rDictionary, const std::string& rFilename){
        std::ofstream file(rFilename, std::ios::out);
        if(file.is_open()){
            file << rDictionary;
            file.close();
        } else {
            QuESo_ERROR << "Warning :: IO::WriteDictionaryToJSON :: Could not create/open file: " << rFilename << ".\n";
        }
    }

    /// @brief Write triangle mesh associated to given condition to STL-file.
    /// @tparam TElementType
    /// @param rCondition
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    template<typename TElementType>
    static void WriteConditionToSTL(const Condition<TElementType>& rCondition,
                                    const std::string& rFilename,
                                    const bool Binary){
        std::ofstream file;
        if(Binary)
            file.open(rFilename, std::ios::out | std::ios::binary);
        else
            file.open(rFilename);

        QuESo_ERROR_IF( !file.is_open() ) << "Could not open file: " << rFilename << ".\n";

        IndexType num_triangles = 0;
        for( const auto& p_segments : rCondition.GetSegments() ){
            num_triangles += p_segments->GetTriangleMesh().NumOfTriangles();
        }

        if(Binary) {
            file << "FileType: Binary                                                                ";
            file.write(reinterpret_cast<const char *>(&static_cast<const int&>(num_triangles)), sizeof(static_cast<const int&>(num_triangles)));

            for( const auto& p_segments : rCondition.GetSegments() ){
                const auto& r_triangle_mesh = p_segments->GetTriangleMesh();
                const IndexType local_num_triangles = r_triangle_mesh.NumOfTriangles();
                for(IndexType triangle_id = 0; triangle_id < local_num_triangles; ++triangle_id) {
                    const auto& p1 = r_triangle_mesh.P1(triangle_id);
                    const auto& p2 = r_triangle_mesh.P2(triangle_id);
                    const auto& p3 = r_triangle_mesh.P3(triangle_id);

                    const auto& normal = r_triangle_mesh.Normal(triangle_id);

                    const float coords[12] = { static_cast<float>(normal[0]), static_cast<float>(normal[1]), static_cast<float>(normal[2]),
                                               static_cast<float>(p1[0]), static_cast<float>(p1[1]), static_cast<float>(p1[2]),
                                               static_cast<float>(p2[0]), static_cast<float>(p2[1]), static_cast<float>(p2[2]),
                                               static_cast<float>(p3[0]), static_cast<float>(p3[1]), static_cast<float>(p3[2]) };

                    for(int i=0; i<12; ++i){
                        file.write(reinterpret_cast<const char *>(&coords[i]), sizeof(coords[i]));
                    }
                    file << "  ";
                }
            }
        }
        else
        {
            file << "solid\n";
            for( const auto& p_segments : rCondition.GetSegments() ){
                const auto& r_triangle_mesh = p_segments->GetTriangleMesh();
                const IndexType local_num_triangles = r_triangle_mesh.NumOfTriangles();

                for(IndexType triangle_id = 0; triangle_id < local_num_triangles; ++triangle_id) {
                    const auto& p1 = r_triangle_mesh.P1(triangle_id);
                    const auto& p2 = r_triangle_mesh.P2(triangle_id);
                    const auto& p3 = r_triangle_mesh.P3(triangle_id);

                    const auto& normal = r_triangle_mesh.Normal(triangle_id);

                    file << "facet normal " << normal[0] << ' ' << normal[1] << ' ' << normal[2] << "\nouter loop\n";
                    file << "vertex " << p1[0] << ' ' << p1[1] << ' ' << p1[2] << ' ' << "\n";
                    file << "vertex " << p2[0] << ' ' << p2[1] << ' ' << p2[2] << ' ' << "\n";
                    file << "vertex " << p3[0] << ' ' << p3[1] << ' ' << p3[2] << ' ' << "\n";
                    file << "endloop\nendfacet\n";
                }
            }
            file << "endsolid"<<std::endl;
        }
        file.close();
    }


    /// @brief Write element container to VTK-file.
    /// @tparam TElementType
    /// @param rBackgroundGrid
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    template<typename TElementType>
    static void WriteElementsToVTK( const BackgroundGrid<TElementType>& rBackgroundGrid,
                                    const std::string& rFilename,
                                    const bool Binary) {
        const SizeType num_elements = rBackgroundGrid.NumberOfActiveElements();

        std::ofstream file;
        if(Binary)
            file.open(rFilename, std::ios::out | std::ios::binary);
        else
            file.open(rFilename);

        QuESo_ERROR_IF( !file.is_open() ) << "Could not open file: " << rFilename << ".\n";

        file << "# vtk DataFile Version 4.1" << std::endl;
        file << "vtk output" << std::endl;
        if(Binary)
            file << "BINARY"<< std::endl;
        else
            file << "ASCII"<< std::endl;

        file << "DATASET UNSTRUCTURED_GRID" << std::endl;
        file << "POINTS " << num_elements*8 << " double" << std::endl;

        const auto begin_el_itr = rBackgroundGrid.ElementsBegin();
        for( IndexType i = 0; i < num_elements; ++i){
            const auto el_itr = (begin_el_itr + i);
            const auto& lower_point = el_itr->GetBoundsXYZ().first;
            const auto& upper_point = el_itr->GetBoundsXYZ().second;

            if( Binary ){
                double rx0 = lower_point[0];
                SwapEnd(rx0);
                double rx1 = upper_point[0];
                SwapEnd(rx1);
                double ry0 = lower_point[1];
                SwapEnd(ry0);
                double ry1 = upper_point[1];
                SwapEnd(ry1);
                double rz0 = lower_point[2];
                SwapEnd(rz0);
                double rz1 = upper_point[2];
                SwapEnd(rz1);

                file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
                file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
                file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

                file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
                file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
                file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

                file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
                file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
                file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

                file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
                file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
                file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

                file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
                file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
                file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

                file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
                file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
                file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

                file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
                file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
                file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

                file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
                file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
                file.write(reinterpret_cast<char*>(&rz1), sizeof(double));
            }
            else {
                file << lower_point[0] << ' ' << lower_point[1] << ' ' << lower_point[2] << std::endl;
                file << upper_point[0] << ' ' << lower_point[1] << ' ' << lower_point[2] << std::endl;
                file << upper_point[0] << ' ' << upper_point[1] << ' ' << lower_point[2] << std::endl;
                file << lower_point[0] << ' ' << upper_point[1] << ' ' << lower_point[2] << std::endl;
                file << lower_point[0] << ' ' << lower_point[1] << ' ' << upper_point[2] << std::endl;
                file << upper_point[0] << ' ' << lower_point[1] << ' ' << upper_point[2] << std::endl;
                file << upper_point[0] << ' ' << upper_point[1] << ' ' << upper_point[2] << std::endl;
                file << lower_point[0] << ' ' << upper_point[1] << ' ' << upper_point[2] << std::endl;
            }
        }
        file << std::endl;
        // Write Cells
        file << "Cells " << num_elements << " " << num_elements*9 << std::endl;
        for( int i = 0; i < static_cast<int>(rBackgroundGrid.NumberOfActiveElements()); ++i){
            if( Binary ){
                int k = 8;
                WriteBinary(file, k);
            for( int j = 0; j < 8; ++j){
                k = 8*i+j;
                WriteBinary(file, k);
            }
            }
            else {
                file << 8 << ' ' << 8*i     << ' ' << 8*i + 1 << ' ' << 8*i + 2 << ' ' << 8*i + 3
                            << ' ' << 8*i + 4 << ' ' << 8*i + 5 << ' ' << 8*i + 6 << ' ' << 8*i + 7 << std::endl;
            }
        }
        file << std::endl;

        file << "CELL_TYPES " << rBackgroundGrid.NumberOfActiveElements() << std::endl;
        for( int i = 0; i < static_cast<int>(rBackgroundGrid.NumberOfActiveElements()); ++i){
            if( Binary ){
                int k = 12;
                WriteBinary(file, k);
            }
            else {
                file << 12 << std::endl;
            }
        }
        file << std::endl;
        file.close();
    }

    /// @brief Write points to VTK. Interface for BackgroundGrid.
    /// @tparam TElementType
    /// @param rBackgroundGrid
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    /// @todo Needs to be refactored.
    template<typename TElementType>
    static void WritePointsToVTK(const BackgroundGrid<TElementType>& rBackgroundGrid,
                                const std::string& rFilename,
                                const bool Binary) {

        const IndexType num_points = rBackgroundGrid.NumberOfIntegrationPoints();
        const IndexType num_elements = num_points;

        std::ofstream file;
        if(Binary)
            file.open(rFilename, std::ios::out | std::ios::binary);
        else
            file.open(rFilename);

        QuESo_ERROR_IF( !file.is_open() ) << "Could not open file: " << rFilename << ".\n";

        file << "# vtk DataFile Version 4.1" << std::endl;
        file << "vtk output" << std::endl;
        if(Binary)
            file << "BINARY"<< std::endl;
        else
            file << "ASCII"<< std::endl;


        file << "DATASET UNSTRUCTURED_GRID" << std::endl;
        file << "POINTS " << num_points << " double" << std::endl;

        const auto el_it_ptr_begin = rBackgroundGrid.ElementsBegin();
        for( IndexType i = 0; i < rBackgroundGrid.NumberOfActiveElements(); ++i){
            const auto el_ptr = (el_it_ptr_begin + i);
            const auto& r_points = el_ptr->GetIntegrationPoints();
            for( const auto& r_point : r_points ){
                auto point_global = el_ptr->PointFromParamToGlobal(r_point.data());

                if( Binary ){
                    WriteBinary(file, point_global[0]);
                    WriteBinary(file, point_global[1]);
                    WriteBinary(file, point_global[2]);
                }
                else {
                    file << point_global[0] << ' ' << point_global[1] << ' ' << point_global[2] << std::endl;
                }
            }

        }
        file << std::endl;

        // Write Cells
        file << "Cells " << num_elements << " " << num_elements*2 << std::endl;
        for( IndexType i = 0; i < num_elements; ++i){
            if( Binary ){
                int k = 1;
                WriteBinary(file, k);
                k = i;
                WriteBinary(file, k);
            }
            else {
                file << 1 << ' ' << i << std::endl;
            }
        }
        file << std::endl;

        file << "CELL_TYPES " << num_elements << std::endl;
        for( IndexType i = 0; i < num_elements; ++i){
            if( Binary ){
                int k = 1;
                WriteBinary(file, k);
            }
            else {
                file << 1 << std::endl;
            }
        }
        file << std::endl;

        file << "POINT_DATA " << num_points << std::endl;
        file << "SCALARS Weights double 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        for( IndexType i = 0; i < rBackgroundGrid.NumberOfActiveElements(); ++i){
            const auto el_ptr = (el_it_ptr_begin + i);
            const auto& points = el_ptr->GetIntegrationPoints();
            for( const auto& point : points ){
                if( Binary ){
                    double rw = point.Weight();
                    WriteBinary(file, rw);
                }
                else {
                    file << point.Weight() << std::endl;
                }
            }
        }
        file << std::endl;
        file.close();
    }

    /// @brief Write points to VTK.
    /// @tparam Type
    /// @param pPoints
    /// @param rFilename
    /// @param Binary If true, file is written in binary format.
    template<typename Type>
    static void WritePointsToVTK(const std::vector<Type>& rPoints,
                                const std::string& rFilename,
                                const bool Binary) {

        const auto begin_points_it_ptr = rPoints.begin();
        const IndexType num_points = rPoints.size();
        const IndexType num_elements = rPoints.size();

        std::ofstream file;
        if(Binary)
            file.open(rFilename, std::ios::out | std::ios::binary);
        else
            file.open(rFilename);

        QuESo_ERROR_IF( !file.is_open() ) << "Could not open file: " << rFilename << ".\n";

        file << "# vtk DataFile Version 4.1" << std::endl;
        file << "vtk output" << std::endl;
        if(Binary)
            file << "BINARY"<< std::endl;
        else
            file << "ASCII"<< std::endl;

        file << "DATASET UNSTRUCTURED_GRID" << std::endl;
        file << "POINTS " << num_points << " double" << std::endl;

        for(IndexType i = 0; i < num_points; ++i){
            auto points_it = (begin_points_it_ptr + i);
            if( Binary ){
                // Make sure to create copy of points. "WriteBinary" will change them.
                auto p_x = (*points_it)[0];
                auto p_y = (*points_it)[1];
                auto p_z = (*points_it)[2];
                WriteBinary(file, p_x);
                WriteBinary(file, p_y);
                WriteBinary(file, p_z);
            }
            else {
                file << (*points_it)[0] << ' ' << (*points_it)[1] << ' ' << (*points_it)[2] << std::endl;
            }
        }
        file << std::endl;

        //Write Cells
        file << "Cells " << num_elements << " " << num_elements*2 << std::endl;
        for( IndexType i = 0; i < num_elements; ++i){
            if( Binary ){
                int k = 1;
                WriteBinary(file, k);
                k = i;
                WriteBinary(file, k);
            }
            else {
                file << 1 << ' ' << i << std::endl;
            }
        }
        file << std::endl;

        file << "CELL_TYPES " << num_elements << std::endl;
        for( IndexType i = 0; i < num_elements; ++i){
            if( Binary ){
                int k = 1;
                WriteBinary(file, k);
            }
            else {
            file << 1 << std::endl;
            }
        }
        file << std::endl;

        file << "POINT_DATA " << num_points << std::endl;
        file << "SCALARS Weights double 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        for(IndexType i = 0; i < num_points; ++i){
            auto points_it = (begin_points_it_ptr + i);

            if( Binary ){
                double rw = points_it->Weight();
                WriteBinary(file, rw);
            }
            else {
                file << points_it->Weight() << std::endl;
            }
        }
        file << std::endl;
        file.close();
    }

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

#endif // IO_UTILTIES_H