#include "utilities/tolerances.h"
#include "containers/triangle_mesh.h"

namespace tibra {

void TriangleMesh::Clean(){
    IndexType size = mTriangles.size();
    IndexType pos = 0;
    while( pos < size ){
        double area = Area(pos);

        if( std::abs(area) < EPS2 ){
            mTriangles.erase( mTriangles.begin() + pos );
            mNormals.erase( mNormals.begin() + pos);
            --size;
        }
        else {
            ++pos;
        }


    }

}
// Protect this Iteration thing in private function.
void TriangleMesh::Refine( IndexType MinNumberOfTriangles ){


    std::vector<IndexType> triangles_ids_to_remove{};

    IndexType original_size = mTriangles.size();
    IndexType size = mTriangles.size();


    double max_area = -1.0;
    double max_edge_length = -1.0;
    for( IndexType pos = 0; pos < size; ++pos){
        max_area = std::max<double>( max_area, Area(pos));
    }

    Reserve(4*size);
    IndexType pos = 0;
    while( pos < size ){
        // Make sure this is a copy!!
        auto vertex_ids = VertexIds(pos);
        const Vector3d p1 = mVertices[vertex_ids[0]];
        const Vector3d p2 = mVertices[vertex_ids[1]];
        const Vector3d p3 = mVertices[vertex_ids[2]];

        const double area = Area(pos);
        if( area > 0.5*max_area ){

            IndexType e1 = AddVertex(CenterEdge(p1, p2));
            IndexType e2 = AddVertex(CenterEdge(p2, p3));
            IndexType e3 = AddVertex(CenterEdge(p3, p1));

            const auto normal = Normal(pos);

            mTriangles.push_back( {vertex_ids[0], e1, e3} );
            mTriangles.push_back( {e1, vertex_ids[1], e2} );
            mTriangles.push_back( {e2, vertex_ids[2], e3} );
            mTriangles.push_back( {e1, e2, e3} );
            mNormals.push_back( normal );
            mNormals.push_back( normal );
            mNormals.push_back( normal );
            mNormals.push_back( normal );
            size += 4;

            mNormals.erase( mNormals.begin() + pos);
            mTriangles.erase(mTriangles.begin() + pos);
            --size;
        } else {
            ++pos;
        }
        if( pos > original_size ){
            Reserve(4*original_size);
            original_size = size;
        }
    }
    if( mTriangles.size() < MinNumberOfTriangles ){
        Refine(MinNumberOfTriangles );
    }
}
}