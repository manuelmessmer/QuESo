// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIANGLE_MESH_INCLUDE_H
#define TRIANGLE_MESH_INCLUDE_H

/// External includes
#include <vector>
#include <array>
#include <cmath>

/// Project includes
#include "geometries/triangle_gauss_legendre_integration_points.h"

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
    typedef std::vector<IntegrationPoint> IntegrationPointVectorType;
    typedef std::unique_ptr<IntegrationPointVectorType> IntegrationPointVectorPtrType;
    ///@}
    ///@name Operations
    ///@{

    /// @brief Area of triangle.
    /// @param TriangleId
    /// @return double.
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

    /// @brief Get integration points in global space.
    /// @param TriangleId
    /// @param Method integration method.
    /// @return IntegrationPointVectorPtrType.
    IntegrationPointVectorPtrType GetIntegrationPointsGlobal( IndexType TriangleId, IndexType Method ) {

        const auto& s_integration_points = GetIntegrationPoints(Method);
        const SizeType point_numbers = s_integration_points.size();

        IntegrationPointVectorPtrType p_global_integration_points
            = std::make_unique<IntegrationPointVectorType>(point_numbers);

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

            (*p_global_integration_points)[i] = IntegrationPoint(xx, yy, zz, s_integration_points[i].GetWeightConst() );
        }

        return std::move(p_global_integration_points);
    }

    const Vector3d& P1(IndexType TriangleId) const {
        return mVertices[mTriangles[TriangleId][0]];
    }

    const Vector3d& P2(IndexType TriangleId) const {
        return mVertices[mTriangles[TriangleId][1]];
    }

    const Vector3d& P3(IndexType TriangleId) const {
        return mVertices[mTriangles[TriangleId][2]];
    }

    void Clear(){
        mVertices.clear();
        mNormals.clear();
        mTriangles.clear();
    }

    void Reserve(IndexType value){
        mVertices.reserve(value);
        mNormals.reserve(value);
        mTriangles.reserve(value);
    }

    void AddVertex(Vector3d NewVertex) {
        return mVertices.push_back(NewVertex);
    }

    void AddTriangle(Vector3i NewTriangle) {
        return mTriangles.push_back(NewTriangle);
    }

    void AddNormal(Vector3d NewNormal) {
        return mNormals.push_back(NewNormal);
    }

    IndexType NumOfTriangles() const{
        return mTriangles.size();
    }

    IndexType NumOfVertices() const{
        return mVertices.size();
    }

    ///@}
private:

    ///@name Private Member Variables
    ///@{
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

    static const std::vector<IntegrationPointVectorType>& AllIntegrationPoints()
    {
        static const std::vector<IntegrationPointVectorType> integration_points =
        {
            TriangleGaussLegendrePoints1::IntegrationPoints(),
            TriangleGaussLegendrePoints2::IntegrationPoints(),
            TriangleGaussLegendrePoints3::IntegrationPoints(),
            TriangleGaussLegendrePoints4::IntegrationPoints()
        };

        return integration_points;
    }

    // Todo: Add checl is method is larger than 3
    static const IntegrationPointVectorType& GetIntegrationPoints( IndexType Method ){

        return AllIntegrationPoints()[Method];
    }

    ///@}
    ///@name Private Member Variables
    ///@{
    std::vector<Vector3d> mVertices;
    std::vector<Vector3i> mTriangles;
    std::vector<Vector3d> mNormals;

    ///@}
};

///@}

#endif