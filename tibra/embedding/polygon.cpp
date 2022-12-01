// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// Project includes
#include "embedding/polygon.h"
#include "utilities/utilities.h"

namespace tibra {

// Function Definitions Polygon
template<IndexType SIZE>
void Polygon<SIZE>::AddVertex(const PointType& rPoint) {
    if( mNumVertices >= SIZE ){
        throw std::runtime_error("Polygon :: AddVertex :: Size of Polygon is exceeded.");
    }
    mVertices[mNumVertices] = rPoint;
    ++mNumVertices;
}

template<IndexType SIZE>
IndexType Polygon<SIZE>::NumVertices() const {
    return mNumVertices;
}

template<IndexType SIZE>
const PointType& Polygon<SIZE>::GetVertex(IndexType i) const {
    if( i >= mNumVertices ){
        throw std::runtime_error("Polygon :: GetVertex :: Size of Polygon is exceeded.");
    }
    return mVertices[i];
}

template<IndexType SIZE>
const PointType& Polygon<SIZE>::operator[] (IndexType i) const {
    return mVertices[i];
}

template<IndexType SIZE>
std::unique_ptr<TriangleMesh> Polygon<SIZE>::pGetTriangleMesh() const {

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

template<IndexType SIZE>
const PointType& Polygon<SIZE>::GetLastVertex() const{
    return mVertices[mNumVertices-1];
}

template<IndexType SIZE>
void Polygon<SIZE>::Clear(){
    std::fill(mVertices.begin(), mVertices.end(), PointType{});
    mNumVertices = 0;
}



template<IndexType SIZE>
std::unique_ptr<typename Polygon<SIZE>::EdgesType> Polygon<SIZE>::pGetEdgesOnPlane(IndexType PlaneIndex, double Position, double PlaneThickness) const {

    auto p_new_edges = std::make_unique<EdgesType>();

    if(mNumVertices < 3){
        throw std::runtime_error("Obacht");
    }

    for( IndexType i = 0 ; i < mNumVertices-1; ++i){
        if( std::abs((mVertices[i][PlaneIndex] - Position)) < PlaneThickness
                && std::abs((mVertices[i+1][PlaneIndex] - Position)) < PlaneThickness ){
            p_new_edges->push_back( {mVertices[i], mVertices[i+1]} );
        }
    }
    if( std::abs((mVertices[mNumVertices-1][PlaneIndex] - Position)) < PlaneThickness
            && std::abs((mVertices[0][PlaneIndex] - Position)) < PlaneThickness ){
        p_new_edges->push_back( {mVertices[mNumVertices-1], mVertices[0]} );
    }

    return std::move(p_new_edges);
}

// Explicit instantiation Polygon
template class Polygon<9>;

} // End namespace tibra