// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef POLYGON_INCLUDE_H
#define POLYGON_INCLUDE_H

//// STL includes
#include <cstddef>
#include <array>
//// Project includes
#include "containers/triangle_mesh.h"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  Polygon
 * @author Manuel Messmer
 * @brief Simple polygon class with static containers.
*/
template<std::size_t SIZE>
class Polygon {

public:
    ///@name Type Definitions
    ///@{

    /// VertexType: Stores PointType and array, that holds information whether or not a points lies on a plane.
    /// Mapping (-x, x, -y, y, -z, z).
    typedef std::pair<PointType, std::array<bool,6> > VertexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    Polygon(const PointType& rNormal) : mNormal(rNormal){}

    ///@}
    ///@name Operations
    ///@{

    ///@brief Adds Point to Polygon. Note that this function sets std::array<bool,6> = false.
    ///       Only use this function, if point on plane information is not required.
    ///@param rPoint PointType
    ///@return Index of added vertex.
    IndexType AddVertex(const PointType& rPoint);

    ///@brief Adds Vertex to Polygon
    ///@param rPoint
    ///@param rIsOnPlanes (-x, +x, -y, y, -z, z) True, if point is on one of respective planes.
    ///@return Index of added vertex.
    IndexType AddVertex(const PointType& rPoint, const std::array<bool, 6>& rIsOnPlanes);

    ///@brief Adds Vertex to Polygon
    ///@param VertexType
    ///@return Index of added vertex.
    IndexType AddVertex(const VertexType& rPoint);

    /// @brief Define point to be on plane.
    /// @param PointIndex
    /// @param PlaneIndex (-x, x, -y, y, -z, z)
    void SetIndexOnPlane(IndexType PointIndex, IndexType PlaneIndex) {
        mVertices[PointIndex].second[PlaneIndex] = true;
    }

    ///@brief Return current number of vertices.
    ///@return IndexType
    IndexType NumVertices() const;

    ///@brief Get i-th vertex.
    ///@param i Index.
    ///@return const VertexType&
    const VertexType& GetVertex(IndexType i) const;

    ///@brief Get i-th vertex.
    ///@param i Index.
    ///@return const VertexType&
    const VertexType& operator[] (IndexType i) const;

    ///@brief Get last vertex.
    ///@return const VertexType&
    const VertexType& GetLastVertex() const;

    ///@brief Triangulates polygon based on centroide and returns triangles. Centroid is computed as mean of all vertices.
    ///@return std::unique_ptr<TriangleMesh> Contains vertices of triangles.
    std::unique_ptr<TriangleMesh> pGetTriangleMesh() const;

    ///@brief Clears vertex container of polygon.
    void Clear();
    ///@}

private:
    ///@name Member variables
    ///@{
    std::array<VertexType, SIZE> mVertices{}; // Keep static array to be fast.
    IndexType mNumVertices = 0;
    PointType mNormal{};
    ///@}
}; // End Polygon class

///@} End TIBRA classes

} // End namespace tibra
#endif // POLYGON_INCLUDE_H