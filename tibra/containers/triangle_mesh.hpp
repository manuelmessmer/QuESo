// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIANGLE_MESH_INCLUDE_H
#define TRIANGLE_MESH_INCLUDE_H

//// STL includes
#include <vector>
#include <array>
#include <cmath>
#include <memory>
#include <iostream>

//// Project includes
#include "define.hpp"
#include "utilities/math_utilities.hpp"
#include "containers/triangle_gauss_legendre_integration_points.hpp"
#include "containers/boundary_integration_point.hpp"

namespace tibra {

///@name TIBRA Classes
///@{
/**
 * @class  TriangleMesh
 * @author Manuel Messmer
 * @brief  Simple implementation of a triangular surface mesh.
 * @todo Refactor AddTriangle() and RemoveTriangle() such that normals are always passed.
 *       Put dot product etc. into math function.
*/
class TriangleMesh
{
public:
    ///@name Type Definitions
    ///@{
    typedef std::vector<IntegrationPoint> IpVectorType;
    typedef Unique<IpVectorType> IpVectorPtrType;
    typedef std::vector<BoundaryIntegrationPoint> BoundaryIpVectorType;
    typedef Unique<BoundaryIpVectorType> BoundaryIpVectorPtrType;
    typedef std::vector<std::vector<std::tuple<IndexType, IndexType, IndexType>>> EdgesOnPlanesVectorType;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Area of triangle.
    /// @details Computation uses cross product.
    /// @note Do not use computation via height of triangle. Not accurate enough.
    /// @param TriangleId
    /// @return double.
    double Area(IndexType TriangleId) const {
        const auto& p1 = this->P1(TriangleId);
        const auto& p2 = this->P2(TriangleId);
        const auto& p3 = this->P3(TriangleId);

        const auto A = p2-p1;
        const auto B = p3-p1;

        return 0.5*Math::Cross(A,B).Norm();
    }

    /// @brief Outward pointing normal.
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

        return (P1+P2+P3) * (1.0/3.0);
    }

    /// @brief Returns aspect ratio.
    /// @param TriangleId
    /// @return double
    double AspectRatio(IndexType TriangleId) const {
        const auto P1 = this->P1(TriangleId);
        const auto P2 = this->P2(TriangleId);
        const auto P3 = this->P3(TriangleId);

        const auto area = this->Area(TriangleId);

        const double a = (P2-P1).Norm(); // length a
        const double b = (P3-P2).Norm(); // length b
        const double c = (P1-P3).Norm(); // length c

        const double max_edge = std::max(std::max(a, b), c);

        const double square_root_3 = 1.73205080756887729;
        return max_edge*(a+b+c)/(4.0*square_root_3 * area);
    }

    /// @brief Get boundary integration points in global space.
    /// @param TriangleId
    /// @param Method integration method.
    /// @return BoundaryIpVectorPtrType.
    BoundaryIpVectorPtrType pGetIPsGlobal( IndexType TriangleId, IndexType Method ) const {

        const auto& s_integration_points = GetIntegrationPoints(Method);
        const SizeType point_numbers = s_integration_points.size();

        auto p_global_integration_points = MakeUnique<BoundaryIpVectorType>(point_numbers);

        const auto& P1 = this->P1(TriangleId);
        const auto& P2 = this->P2(TriangleId);
        const auto& P3 = this->P3(TriangleId);

        for( int i = 0; i < point_numbers; ++i){
            const double xx  = ShapeFunctionValue( 0, s_integration_points[i] ) * P1[0] +
                               ShapeFunctionValue( 1, s_integration_points[i] ) * P2[0] +
                               ShapeFunctionValue( 2, s_integration_points[i] ) * P3[0] ;

            const double yy = ShapeFunctionValue( 0, s_integration_points[i] ) * P1[1] +
                              ShapeFunctionValue( 1, s_integration_points[i] ) * P2[1] +
                              ShapeFunctionValue( 2, s_integration_points[i] ) * P3[1] ;

            const double zz = ShapeFunctionValue( 0, s_integration_points[i] ) * P1[2] +
                              ShapeFunctionValue( 1, s_integration_points[i] ) * P2[2] +
                              ShapeFunctionValue( 2, s_integration_points[i] ) * P3[2] ;

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
        mEdgesOnPlanes.clear();
    }

    ///@brief Reserve all containers for normal vertices and triangles.
    ///       Note, call ReserveEdgesOnPlane() to reserve edges containers.
    ///@param Size
    void Reserve(IndexType Size){
        mVertices.reserve(Size);
        mNormals.reserve(Size);
        mTriangles.reserve(Size);
        /// Resize to avoid segmenation faults, if ReserveEdgesOnPlane() is not called.
        if( mEdgesOnPlanes.size() != 6 ){
            mEdgesOnPlanes.resize(6);
        }
    }

    /// @brief Reserve containers for edges on plane.
    /// @param Size
    void ReserveEdgesOnPlane(IndexType Size){
        if( mEdgesOnPlanes.size() != 6 ){
            mEdgesOnPlanes.resize(6);
        }
        for( auto& edges : mEdgesOnPlanes ){
            edges.reserve(Size);
        }
    }

    ///@brief Add vertex to mesh.
    ///@param NewVertex
    IndexType AddVertex(const Vector3d& NewVertex) {
        mVertices.push_back(NewVertex);
        return mVertices.size()-1;
    }

    ///@brief Add triangle to mesh.
    ///@param NewTriangle
    void AddTriangle(const Vector3i& NewTriangle) {
        mTriangles.push_back(NewTriangle);
    }

    /// @brief Remove triangle by index.
    /// @param Index
    void RemoveTriangle(IndexType Index ){
        mTriangles.erase( mTriangles.begin() + Index );
    }

    /// @brief Remove normal by index.
    /// @param Index
    void RemoveNormal(IndexType Index ){
        mNormals.erase( mNormals.begin() + Index );
    }

    ///@brief Add normal to mesh.
    ///@param NewNormal
    void AddNormal(const Vector3d& NewNormal) {
        mNormals.push_back(NewNormal);
    }

    /// @brief Add edges on plane to triangle mesh. Note, edges on plane are only required for clipped meshes.
    /// @param PlaneIndex Index of plane [-x, x, -y, y, -z, z].
    /// @param V1 Index Vertex 1.
    /// @param V2 Index Vertex 2.
    /// @param Normal Normal of triangle attached to edge.
    void AddEdgeOnPlane(IndexType PlaneIndex, IndexType V1, IndexType V2, IndexType Normal){
        mEdgesOnPlanes[PlaneIndex].push_back( std::make_tuple(V1, V2, Normal) );
    }

    /// @brief Returns all edges on planes.
    /// @return EdgesOnPlanesVectorType&
    const EdgesOnPlanesVectorType& GetEdgesOnPlanes() const {
        return mEdgesOnPlanes;
    }
    ///@brief Get number of triangles in mesh.
    IndexType NumOfTriangles() const{
        return mTriangles.size();
    }

    ///@brief Get number of vertices in mesh.
    IndexType NumOfVertices() const{
        return mVertices.size();
    }

    ///@brief Get vertices from mesh. (const version)
    ///@return const std::vector<Vector3d>&
    const std::vector<Vector3d>& GetVertices() const {
        return mVertices;
    }

    /// @brief Get vertices from mesh. (non-const version)
    /// @return std::vector<Vector3d>&
    std::vector<Vector3d>& GetVertices() {
        return mVertices;
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

    ///@name Private Member Variables
    ///@{

    ///@brief Returns ShapeFunctionValue
    ///@param ShapeFunctionIndex
    ///@param rPoint
    ///@return double
    static double ShapeFunctionValue( IndexType ShapeFunctionIndex, const Vector3d& rPoint ) {
        switch( ShapeFunctionIndex )
        {
        case 0:
            return( 1.0 -rPoint[0] - rPoint[1] );
        case 1:
            return( rPoint[0] );
        case 2:
            return( rPoint[1] );
        default:
            TIBRA_ERROR("TriangleMesh::ShapeFunctionValue") << "Wrong Index of Shape Function.\n";
            break;
        }

        return 0;
    }

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
        TIBRA_ERROR_IF("TriangleMesh::GetIntegrationPoints", Method > 3) << "IntegrationPoint Index exceeds default.\n";

        return AllIntegrationPoints()[Method];
    }

    ///@}
    ///@name Private Member Variables
    ///@{
    std::vector<Vector3d> mVertices;
    std::vector<Vector3i> mTriangles;
    std::vector<Vector3d> mNormals;
    EdgesOnPlanesVectorType mEdgesOnPlanes{};

    ///@}

}; // End of class TriangleMesh
///@} // End TIBRA classes

} // End namespace tibra

#endif // TRIANGLE_MESH_INCLUDE_H