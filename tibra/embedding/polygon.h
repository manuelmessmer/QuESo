// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef POLYGON_INCLUDE_H
#define POLYGON_INCLUDE_H

/// External libraries
#include <cstddef>
#include <array>

///Project includes
#include "containers/triangle_mesh.h"

///@name TIBRA Classes
///@{

/**
 * @class  Polygon
 * @author Manuel Messmer
 * @brief Simple polygon class with static containers.
*/
template<std::size_t DIM, std::size_t SIZE>
class Polygon {

public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType;
    typedef std::array<double, DIM> PointType;
    typedef std::vector<std::array<PointType, 2>> EdgesType;

    Polygon(const PointType& rNormal) : mNormal(rNormal){}

    ///@}
    ///@name Operations
    ///@{

    ///@brief Adds Vertex to Polygon
    ///@param rPoint
    void AddVertex(const PointType& rPoint);

    ///@brief Return current number of vertices.
    ///@return IndexType
    IndexType NumVertices() const;

    ///@brief Get i-th vertex.
    ///@param i Index.
    ///@return const PointType&
    const PointType& GetVertex(IndexType i) const;

    ///@brief Get i-th vertex.
    ///@param i Index.
    ///@return const PointType&
    const PointType& operator[] (IndexType i) const;

    ///@brief Get last vertex.
    const PointType& GetLastVertex() const;

    ///@brief Triangulates polygon centroid coodinate and returns triangles. Centroid is computed as mean of all vertices.
    ///@param rNormal Provide normal to avoid recomputation.
    ///@return std::unique_ptr<TriangleMesh> Contains vertices of triangles.
    std::unique_ptr<TriangleMesh> pGetTriangleMesh() const;

    ///@brief GetBoundaryEdges
    std::unique_ptr<EdgesType> pGetEdgesOnPlane(IndexType PlaneIndex, double Position, double PlaneThickness) const;

    ///@brief Clears vertex container of polygon.
    void Clear();
    ///@}

private:
    ///@name Member variables
    ///@{
    std::array<PointType, SIZE> mVertices{}; // Keep static array to be fast.
    IndexType mNumVertices = 0;
    PointType mNormal{};
    ///@}
};

///@}

#endif