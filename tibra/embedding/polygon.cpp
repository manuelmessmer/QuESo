// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// Project includes
#include "embedding/polygon.h"
#include "utilities/utilities.h"

namespace tibra {

// Function Definitions Polygon
template<IndexType SIZE>
IndexType Polygon<SIZE>::AddVertex(const PointType& rPoint ) {
    if( mNumVertices >= SIZE ){
        throw std::runtime_error("Polygon :: AddVertex :: Size of Polygon is exceeded.");
    }
    mVertices[mNumVertices].first = rPoint;
    return mNumVertices++;
}

template<IndexType SIZE>
IndexType Polygon<SIZE>::AddVertex(const PointType& rPoint, const std::array<bool, 6>& rStatus) {
    if( mNumVertices >= SIZE ){
        throw std::runtime_error("Polygon :: AddVertex :: Size of Polygon is exceeded.");
    }
    mVertices[mNumVertices].first = rPoint;
    mVertices[mNumVertices].second = rStatus;
    return mNumVertices++;
}

template<IndexType SIZE>
IndexType Polygon<SIZE>::AddVertex(const VertexType& rPoint) {
    if( mNumVertices >= SIZE ){
        throw std::runtime_error("Polygon :: AddVertex :: Size of Polygon is exceeded.");
    }
    mVertices[mNumVertices] = rPoint;
    return mNumVertices++;
}

template<IndexType SIZE>
IndexType Polygon<SIZE>::NumVertices() const {
    return mNumVertices;
}

template<IndexType SIZE>
const typename Polygon<SIZE>::VertexType& Polygon<SIZE>::GetVertex(IndexType i) const {
    if( i >= mNumVertices ){
        throw std::runtime_error("Polygon :: GetVertex :: Size of Polygon is exceeded.");
    }
    return mVertices[i];
}

template<IndexType SIZE>
const typename Polygon<SIZE>::VertexType& Polygon<SIZE>::operator[] (IndexType i) const {
    if( i >= mNumVertices ){
        throw std::runtime_error("Polygon :: GetVertex :: Size of Polygon is exceeded.");
    }
    return mVertices[i];
}

template<IndexType SIZE>
const typename Polygon<SIZE>::VertexType& Polygon<SIZE>::GetLastVertex() const{
    return mVertices[mNumVertices-1];
}

template<IndexType SIZE>
void Polygon<SIZE>::Clear(){
    const std::array<bool,6> points_on_plane = {false, false, false, false, false, false};
    std::fill(mVertices.begin(), mVertices.end(), std::make_pair(PointType{}, points_on_plane) );
    mNumVertices = 0;
}

template<IndexType SIZE>
std::unique_ptr<TriangleMesh> Polygon<SIZE>::pGetTriangleMesh() const {
    auto p_new_mesh = std::make_unique<TriangleMesh>();
    p_new_mesh->Reserve(mNumVertices);

    if(mNumVertices < 3){
        return std::move(p_new_mesh);
    }

    if( mNumVertices == 3 ){
        p_new_mesh->AddVertex( mVertices[0].first );
        p_new_mesh->AddVertex( mVertices[1].first );
        p_new_mesh->AddVertex( mVertices[2].first );

        p_new_mesh->AddTriangle( {0, 1, 2} );
        p_new_mesh->AddNormal( mNormal );

        // Add edges, that are located on a plane, to the mesh.
        // Planes: (-x, +x, -y, y, -z, z)
        for( IndexType plane_index = 0; plane_index < 6; ++plane_index ){
            const bool v1_on_plane = mVertices[0].second[plane_index];
            const bool v2_on_plane = mVertices[1].second[plane_index];
            const bool v3_on_plane = mVertices[2].second[plane_index];
            if( (v1_on_plane + v2_on_plane + v3_on_plane) == 3 ) {
                throw std::runtime_error("Polygon :: pGetTriangleMesh :: All vertices are set on plane.");
            }
            if( v1_on_plane && v2_on_plane ){
                p_new_mesh->AddEdgeOnPlane(plane_index, 0, 1, 0);
            }
            else if( v2_on_plane && v3_on_plane ){
                p_new_mesh->AddEdgeOnPlane(plane_index, 1, 2, 0);
            }
            else if( v3_on_plane && v1_on_plane ){
                p_new_mesh->AddEdgeOnPlane(plane_index, 2, 0, 0);
            }
        }

        return std::move(p_new_mesh);
    }

    // Compute mean of vertices
    PointType centroid = {0.0, 0.0, 0.0};
    for( IndexType i = 0 ; i < mNumVertices; ++i){
        centroid += mVertices[i].first;
    }
    centroid /= mNumVertices;

    IndexType vertex_count = 0;
    IndexType triangle_count = 0;
    for( IndexType i = 0 ; i < mNumVertices-1; ++i){
        p_new_mesh->AddVertex( mVertices[i].first );
        p_new_mesh->AddVertex( mVertices[i+1].first );
        p_new_mesh->AddVertex( centroid );

        p_new_mesh->AddTriangle( {vertex_count+0, vertex_count+1, vertex_count+2} );
        p_new_mesh->AddNormal( mNormal );

        for( IndexType plane_index = 0; plane_index < 6; ++plane_index ){
            const bool v1_on_plane = mVertices[i].second[plane_index];
            const bool v2_on_plane = mVertices[i+1].second[plane_index];
            if( v1_on_plane && v2_on_plane ){
                p_new_mesh->AddEdgeOnPlane(plane_index, vertex_count+0, vertex_count+1, triangle_count);
            }
        }
        ++triangle_count;
        vertex_count += 3;
    }

    p_new_mesh->AddVertex( mVertices[mNumVertices-1].first );
    p_new_mesh->AddVertex( mVertices[0].first );
    p_new_mesh->AddVertex( centroid );

    p_new_mesh->AddTriangle( {vertex_count+0, vertex_count+1, vertex_count+2 } );
    p_new_mesh->AddNormal( mNormal );

    for( IndexType plane_index = 0; plane_index < 6; ++plane_index ){
        const bool v1_on_plane = mVertices[mNumVertices-1].second[plane_index];
        const bool v2_on_plane = mVertices[0].second[plane_index];
        if( v1_on_plane && v2_on_plane ){
            p_new_mesh->AddEdgeOnPlane(plane_index, vertex_count+0, vertex_count+1, triangle_count);
        }
    }

    return std::move(p_new_mesh);
}


// Explicit instantiation Polygon
template class Polygon<9>;
template class Polygon<4>;

} // End namespace tibra