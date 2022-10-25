// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// Project includes
#include "embedding/polygon.h"
#include "utilities/utilities.h"

typedef std::size_t IndexType;

// Function Definitions Polygon
template<std::size_t DIM, IndexType SIZE>
void Polygon<DIM, SIZE>::AddVertex(const PointType& rPoint) {
    if( mNumVertices >= SIZE ){
        throw std::runtime_error("Polygon :: AddVertex :: Size of Polygon is exceeded.");
    }
    mVertices[mNumVertices] = rPoint;
    ++mNumVertices;
}

template<std::size_t DIM, IndexType SIZE>
IndexType Polygon<DIM, SIZE>::NumVertices() const {
    return mNumVertices;
}

template<std::size_t DIM, IndexType SIZE>
const typename Polygon<DIM, SIZE>::PointType& Polygon<DIM, SIZE>::GetVertex(IndexType i) const {
    if( i >= mNumVertices ){
        throw std::runtime_error("Polygon :: GetVertex :: Size of Polygon is exceeded.");
    }
    return mVertices[i];
}

template<std::size_t DIM, IndexType SIZE>
const typename Polygon<DIM, SIZE>::PointType& Polygon<DIM, SIZE>::operator[] (IndexType i) const {
    return mVertices[i];
}

template<std::size_t DIM, IndexType SIZE>
std::unique_ptr<TriangleMesh> Polygon<DIM, SIZE>::pGetTriangleMesh() const {

    auto p_new_mesh = std::make_unique<TriangleMesh>();
    p_new_mesh->Reserve(mNumVertices);

    if(mNumVertices < 3){
        throw std::runtime_error("Clipper :: pGetTriangleMesh :: mNumVertices < 3.");
    }

    if( mNumVertices == 3 ){

        p_new_mesh->AddVertex( mVertices[0] );
        p_new_mesh->AddVertex( mVertices[1] );
        p_new_mesh->AddVertex( mVertices[2] );

        p_new_mesh->AddTriangle( {0, 1, 2} );
        p_new_mesh->AddNormal( mNormal );

        return std::move(p_new_mesh);
    }

    // Compute mean of vertices
    PointType centroid = {0.0, 0.0, 0.0};
    const double inv_num_vertices = 1.0/mNumVertices;

    for( IndexType i = 0 ; i < mNumVertices; ++i){
        centroid[0] += inv_num_vertices*mVertices[i][0];
        centroid[1] += inv_num_vertices*mVertices[i][1];
        centroid[2] += inv_num_vertices*mVertices[i][2];
    }

    IndexType vertex_count = 0;
    for( IndexType i = 0 ; i < mNumVertices-1; ++i){
        p_new_mesh->AddVertex( mVertices[i] );
        p_new_mesh->AddVertex( mVertices[i+1] );
        p_new_mesh->AddVertex( centroid );

        p_new_mesh->AddTriangle( {vertex_count+0, vertex_count+1, vertex_count+2} );
        p_new_mesh->AddNormal( mNormal );
        vertex_count += 3;
    }

    p_new_mesh->AddVertex( mVertices[mNumVertices-1] );
    p_new_mesh->AddVertex( mVertices[0] );
    p_new_mesh->AddVertex( centroid );

    p_new_mesh->AddTriangle( {vertex_count+0, vertex_count+1, vertex_count+2 } );
    p_new_mesh->AddNormal( mNormal );

    return std::move(p_new_mesh);
}

template<std::size_t DIM, IndexType SIZE>
const typename Polygon<DIM, SIZE>::PointType& Polygon<DIM, SIZE>::GetLastVertex() const{
    return mVertices[mNumVertices-1];
}

template<std::size_t DIM, IndexType SIZE>
void Polygon<DIM, SIZE>::Clear(){
    std::fill(mVertices.begin(), mVertices.end(), PointType{});
    mNumVertices = 0;
}

template<std::size_t DIM, IndexType SIZE>
std::unique_ptr<BoundaryEdges<0,1>> Polygon<DIM, SIZE>::pGetBoundaryEdges(IndexType PlaneIndex, double Position, double PlaneThickness) const {

    auto p_new_edges = std::make_unique<BoundaryEdges<0,1>>();

    if(mNumVertices < 3){
        throw std::runtime_error("Obacht");
    }

    if( mNumVertices == 3 ){
        return nullptr;
    }

    for( IndexType i = 0 ; i < mNumVertices-1; ++i){
        if( std::abs((mVertices[i][PlaneIndex] - Position))
                && std::abs((mVertices[i+1][PlaneIndex] - Position)) ){
            p_new_edges->AddEdge(mVertices[i], mVertices[i+1], mNormal);
        }
    }
    if( std::abs((mVertices[mNumVertices-1][PlaneIndex] - Position))
            && std::abs((mVertices[0][PlaneIndex] - Position)) ){
        p_new_edges->AddEdge(mVertices[mNumVertices-1], mVertices[0], mNormal );
    }

    return std::move(p_new_edges);;
}

// Explicit instantiation Polygon
template class Polygon<3, 9>;