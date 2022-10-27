// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIANGLE_MESH_INCLUDE_H
#define TRIANGLE_MESH_INCLUDE_H

/// External includes
#include <vector>
#include <array>
#include <cmath>
#include <memory>
#include <iostream>
#include <map>

/// Project includes
#include "geometries/triangle_gauss_legendre_integration_points.h"
#include "geometries/boundary_integration_point.h"

///@name TIBRA Classes
///@{

///
/**
 * @class  TriangleMesh
 * @author Manuel Messmer
 * @brief  Simple implementation of a triangular surface mesh.
*/
class TriangleMesh
{
public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef std::array<double,3> Vector3d;
    typedef std::array<IndexType,3> Vector3i;
    typedef std::vector<IntegrationPoint> IpVectorType;
    typedef std::unique_ptr<IpVectorType> IpVectorPtrType;
    typedef std::vector<BoundaryIntegrationPoint> BoundaryIpVectorType;
    typedef std::unique_ptr<BoundaryIpVectorType> BoundaryIpVectorPtrType;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Area of triangle.
    /// @param TriangleId
    /// @return double.
    /// @todo Consider using more efficient std::pow alternative.
    double Area(IndexType TriangleId) const {
        const auto P1 = this->P1(TriangleId);
        const auto P2 = this->P2(TriangleId);
        const auto P3 = this->P3(TriangleId);

        const double a = std::sqrt( std::pow(P1[0] - P2[0], 2) + std::pow(P1[1] - P2[1], 2) + std::pow(P1[2] - P2[2], 2));
        const double b = std::sqrt( std::pow(P2[0] - P3[0], 2) + std::pow(P2[1] - P3[1], 2) + std::pow(P2[2] - P3[2], 2));
        const double c = std::sqrt( std::pow(P3[0] - P1[0], 2) + std::pow(P3[1] - P1[1], 2) + std::pow(P3[2] - P1[2], 2));

        const double s = (a+b+c) / 2.0;

        return std::sqrt(s*(s-a)*(s-b)*(s-c));
    }

    /// @brief Outward point normal.
    /// @param TriangleId
    /// @return Vector3d.
    const Vector3d& Normal(IndexType TriangleId) const{
        return mNormals[TriangleId];
    }

    /// @brief Center of triangles in global coordinates.
    /// @param TriangleId
    /// @return Vector3d.
    Vector3d Center(IndexType TriangleId) const {
        const auto P1 = this->P1(TriangleId);
        const auto P2 = this->P2(TriangleId);
        const auto P3 = this->P3(TriangleId);

        Vector3d tmp_point = {0.0, 0.0, 0.0};
        tmp_point[0] = 1.0/3.0 * (P1[0] + P2[0] + P3[0]);
        tmp_point[1] = 1.0/3.0 * (P1[1] + P2[1] + P3[1]);
        tmp_point[2] = 1.0/3.0 * (P1[2] + P2[2] + P3[2]);

        return tmp_point;
    }

    /// @brief Get boundary integration points in global space.
    /// @param TriangleId
    /// @param Method integration method.
    /// @return IntegrationPointVectorPtrType.
    BoundaryIpVectorPtrType GetIPsGlobal( IndexType TriangleId, IndexType Method ) const {

        const auto& s_integration_points = GetIntegrationPoints(Method);
        const SizeType point_numbers = s_integration_points.size();

        auto p_global_integration_points = std::make_unique<BoundaryIpVectorType>(point_numbers);

        const auto P1 = this->P1(TriangleId);
        const auto P2 = this->P2(TriangleId);
        const auto P3 = this->P3(TriangleId);

        for( int i = 0; i < point_numbers; ++i){
            const double xx  = this->ShapeFunctionValue( 0, s_integration_points[i] ) * P1[0] +
                               this->ShapeFunctionValue( 1, s_integration_points[i] ) * P2[0] +
                               this->ShapeFunctionValue( 2, s_integration_points[i] ) * P3[0] ;

            const double yy = this->ShapeFunctionValue( 0, s_integration_points[i] ) * P1[1] +
                              this->ShapeFunctionValue( 1, s_integration_points[i] ) * P2[1] +
                              this->ShapeFunctionValue( 2, s_integration_points[i] ) * P3[1] ;

            const double zz = this->ShapeFunctionValue( 0, s_integration_points[i] ) * P1[2] +
                              this->ShapeFunctionValue( 1, s_integration_points[i] ) * P2[2] +
                              this->ShapeFunctionValue( 2, s_integration_points[i] ) * P3[2] ;

            // Normalize weights to 1 by multiplying by 2.
            const double weight = 2.0*s_integration_points[i].GetWeight()*Area(TriangleId);
            (*p_global_integration_points)[i] = BoundaryIntegrationPoint(xx, yy, zz, weight, Normal(TriangleId) );
        }

        return std::move(p_global_integration_points);
    }

    ///@brief Get triangle vertex 1
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d& P1(IndexType TriangleId) const {
        return mVertices[mTriangles[TriangleId][0]];
    }

    ///@brief Get triangle vertex 2
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d& P2(IndexType TriangleId) const {
        return mVertices[mTriangles[TriangleId][1]];
    }

    ///@brief Get triangle vertex 3
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d& P3(IndexType TriangleId) const {
        return mVertices[mTriangles[TriangleId][2]];
    }

    ///@brief Get triangle vertex 3
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3i& VertexIds(IndexType TriangleId) const {
        return mTriangles[TriangleId];
    }

    ///@brief Clear all containers.
    void Clear(){
        mVertices.clear();
        mNormals.clear();
        mTriangles.clear();
    }

    ///@brief Reserve all containers.
    ///@param Size
    void Reserve(IndexType Size){
        mVertices.reserve(Size);
        mNormals.reserve(Size);
        mTriangles.reserve(Size);
    }

    ///@brief Add vertex to mesh.
    ///@param NewVertex
    void AddVertex(const Vector3d& NewVertex) {
        return mVertices.push_back(NewVertex);
    }

    ///@brief Add triangle to mesh.
    ///@param NewTriangle
    void AddTriangle(const Vector3i& NewTriangle) {
        return mTriangles.push_back(NewTriangle);
    }

    ///@brief Add normal to mesh.
    ///@param NewNormal
    void AddNormal(const Vector3d& NewNormal) {
        return mNormals.push_back(NewNormal);
    }

    ///@brief Get number of triangles in mesh.
    IndexType NumOfTriangles() const{
        return mTriangles.size();
    }

    ///@brief Get number of vertices in mesh.
    IndexType NumOfVertices() const{
        return mVertices.size();
    }

    ///@brief Get vertices from mesh.
    const std::vector<Vector3d>& GetVertices() const {
        return mVertices;
    }

    // Vector3d GetPointGlobalSpace(IndexType TriangleId, const Vector3d& rPoint){
    //     const auto P1 = this->P1(TriangleId);
    //     const auto P2 = this->P2(TriangleId);
    //     const auto P3 = this->P3(TriangleId);

    //     const double xx = this->ShapeFunctionValue( 0, rPoint ) * P1[0] +
    //                       this->ShapeFunctionValue( 1, rPoint ) * P2[0] +
    //                       this->ShapeFunctionValue( 2, rPoint ) * P3[0] ;

    //     const double yy = this->ShapeFunctionValue( 0, rPoint ) * P1[1] +
    //                       this->ShapeFunctionValue( 1, rPoint ) * P2[1] +
    //                       this->ShapeFunctionValue( 2, rPoint ) * P3[1] ;

    //     const double zz = this->ShapeFunctionValue( 0, rPoint ) * P1[2] +
    //                       this->ShapeFunctionValue( 1, rPoint ) * P2[2] +
    //                       this->ShapeFunctionValue( 2, rPoint ) * P3[2] ;

    //     return Vector3d{xx, yy, zz};
    // }

    ///@brief Copy from other mesh.
    ///@param rTriangleMesh Triangle mesh to copy from.
    void Append( const TriangleMesh& rTriangleMesh ){
        std::vector<IndexType> indices{};
        indices.reserve(rTriangleMesh.NumOfTriangles());
        for( IndexType i = 0; i < rTriangleMesh.NumOfTriangles(); ++i){
            indices.push_back(i);
        }
        Append(indices, rTriangleMesh);
    }

    ///@brief Partial copy from other mesh.
    ///@param rTriangleIndices Indices of triangles to be copied.
    ///@param rTriangleMesh Triangle mesh to copy from.
    void Append( const std::vector<IndexType>& rTriangleIndices, const TriangleMesh& rTriangleMesh ){

        IndexType vertex_count = NumOfVertices();
        std::map<IndexType, IndexType> index_map{};

        //const IndexType initial_num_vertices = NumOfVertices();

        for( auto triangle : rTriangleIndices){
            const auto& tmp_indices = rTriangleMesh.VertexIds(triangle);
            std::array<IndexType,3> new_triangle{};
            IndexType ii = 0;
            for( auto index : tmp_indices ){
                // Insert index into index_map if map does not contain index.
                auto ret = index_map.insert( std::pair<IndexType,IndexType>(index, vertex_count) );
                if (ret.second==true) {
                    new_triangle[ii] = vertex_count;
                    vertex_count++;
                } else {
                    new_triangle[ii] = index_map[index];
                }
                ii++;
            }
            // Copy triangles and normals.
            AddTriangle(new_triangle);
            AddNormal( rTriangleMesh.Normal(triangle) );
        }

        mVertices.resize(vertex_count);
        const auto& new_vertices = rTriangleMesh.GetVertices();

        // Copy vertices.
        for( auto index : index_map ){
            mVertices[ index.second ] = new_vertices[ index.first ];
        }

        Check();
    }


    ///@brief Basic check of this TriangleMesh instance.
    bool Check() const{
        // Check if mTriangles and mNormals are of the same size.
        if( mTriangles.size() != mNormals.size() ){
            std::cerr << "TriangleMesh :: Number of Triangles and Normals in mesh do not match.\n";
            return false;
        }
        // Check if all vertex ids exist.
        for( int i = 0; i < mTriangles.size(); ++i ){
            for(int j = 0; j < 3; ++j){
                if( mTriangles[i][j] >= mVertices.size() ){
                    std::cerr << "TriangleMesh :: Triangle/Vertex mismatch.\n";
                    return false;
                }
            }
        }
        return true;
    }
    ///@}
private:

    ///@brief Returns ShapeFunctionValue
    ///@param ShapeFunctionIndex
    ///@param rPoint
    ///@return double
    double ShapeFunctionValue( IndexType ShapeFunctionIndex, const Vector3d& rPoint ) const {
        switch( ShapeFunctionIndex )
        {
        case 0:
            return( 1.0 -rPoint[0] - rPoint[1] );
        case 1:
            return( rPoint[0] );
        case 2:
            return( rPoint[1] );
        default:
            throw std::invalid_argument(" Triangle3D3N Wrong Index of Shape Function! ");
            break;
        }

        return 0;
    }
    ///@name Private Member Variables
    ///@{


    ///@brief Factory function for triangle Gauss Legendre points.
    static const std::vector<IpVectorType>& AllIntegrationPoints()
    {
        static const std::vector<IpVectorType> integration_points =
        {
            TriangleGaussLegendrePoints1::IntegrationPoints(),
            TriangleGaussLegendrePoints2::IntegrationPoints(),
            TriangleGaussLegendrePoints3::IntegrationPoints(),
            TriangleGaussLegendrePoints4::IntegrationPoints()
        };

        return integration_points;
    }

    ///@brief Get triangle Gauss Legendre points by Method - options (0,1,2,3)
    static const IpVectorType& GetIntegrationPoints( IndexType Method ){
        if( Method > 3){
            throw std::runtime_error("TriangleMesh::GetIntegrationPoints IntegrationPoint Index exceeds default.");
        }
        return AllIntegrationPoints()[Method];
    }

    ///@}
    ///@name Private Member Variables
    ///@{
    std::vector<Vector3d> mVertices;
    std::vector<Vector3i> mTriangles;
    std::vector<Vector3d> mNormals;

    ///@}

}; // End of class TriangleMesh

///@}

#endif