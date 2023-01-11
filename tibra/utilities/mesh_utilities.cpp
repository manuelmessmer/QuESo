// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <map>
#include <numeric>
/// Project includes
#include "utilities/mesh_utilities.h"
#include "utilities/math_utilities.hpp"

namespace tibra {

typedef MeshUtilities::TriangleMeshPtrType TriangleMeshPtrType;

void MeshUtilities::Refine(TriangleMesh& rTriangleMesh, IndexType MinNumberOfTriangles){
    IndexType original_size = rTriangleMesh.NumOfTriangles();
    IndexType size = original_size;

    double max_area = -1.0;
    for( IndexType pos = 0; pos < size; ++pos){
        max_area = std::max<double>( max_area, rTriangleMesh.Area(pos));
    }

    const auto& r_vertices = rTriangleMesh.GetVertices();
    rTriangleMesh.Reserve(4*size);
    IndexType pos = 0;
    while( pos < size ){
        // Make sure this is a copy!!
        auto vertex_ids = rTriangleMesh.VertexIds(pos);
        const Vector3d p1 = r_vertices[vertex_ids[0]];
        const Vector3d p2 = r_vertices[vertex_ids[1]];
        const Vector3d p3 = r_vertices[vertex_ids[2]];

        const double area = rTriangleMesh.Area(pos);
        if( area > 0.5*max_area ){
            IndexType e1 = rTriangleMesh.AddVertex( (p1 + p2)*0.5 );
            IndexType e2 = rTriangleMesh.AddVertex( (p2 + p3)*0.5 );
            IndexType e3 = rTriangleMesh.AddVertex( (p3 + p1)*0.5 );

            const auto normal = rTriangleMesh.Normal(pos);
            rTriangleMesh.AddTriangle( {vertex_ids[0], e1, e3} );
            rTriangleMesh.AddTriangle( {e1, vertex_ids[1], e2} );
            rTriangleMesh.AddTriangle( {e2, vertex_ids[2], e3} );
            rTriangleMesh.AddTriangle( {e1, e2, e3} );
            rTriangleMesh.AddNormal( normal );
            rTriangleMesh.AddNormal( normal );
            rTriangleMesh.AddNormal( normal );
            rTriangleMesh.AddNormal( normal );
            size += 4;

            rTriangleMesh.RemoveTriangle(pos);
            rTriangleMesh.RemoveNormal(pos);
            --size;
        } else {
            ++pos;
        }
        if( pos > original_size ){
            rTriangleMesh.Reserve(4*original_size);
            original_size = size;
        }
    }
    if( rTriangleMesh.NumOfTriangles() < MinNumberOfTriangles ){
        Refine(rTriangleMesh, MinNumberOfTriangles );
    }
}

void MeshUtilities::Append(TriangleMesh& rTriangleMesh, const TriangleMesh& rNewMesh){
    std::vector<IndexType> indices(rNewMesh.NumOfTriangles());
    /// Fill vector with number of increasing values: 0,1,2..
    std::iota( indices.begin(), indices.end(), 0 );
    Append(rTriangleMesh, rNewMesh, indices);
}

void MeshUtilities::Append(TriangleMesh& rTriangleMesh, const TriangleMesh& rNewMesh, const std::vector<IndexType>& rIndices){

    const IndexType initial_number_triangles = rTriangleMesh.NumOfTriangles();
    const IndexType initial_number_vertices = rTriangleMesh.NumOfVertices();
    IndexType vertex_count = initial_number_vertices;
    std::map<IndexType, IndexType> index_map_vertices{};

    for( auto triangle : rIndices){
        const auto& tmp_indices = rNewMesh.VertexIds(triangle);
        Vector3i new_triangle{};
        IndexType ii = 0;
        for( auto index : tmp_indices ){
            // Insert index into index_map_vertices if map does not contain index.
            auto ret = index_map_vertices.insert( std::pair<IndexType,IndexType>(index, vertex_count) );
            if (ret.second==true) {
                new_triangle[ii] = vertex_count;
                vertex_count++;
            } else {
                new_triangle[ii] = index_map_vertices[index];
            }
            ++ii;
        }
        // Copy triangles and normals.
        rTriangleMesh.AddTriangle(new_triangle);
        rTriangleMesh.AddNormal( rNewMesh.Normal(triangle) );
    }

    auto& r_vertices = rTriangleMesh.GetVertices();
    r_vertices.resize(vertex_count);
    const auto& r_new_vertices = rNewMesh.GetVertices();

    // Copy vertices.
    for( auto index : index_map_vertices ){
        r_vertices[ index.second ] = r_new_vertices[ index.first ];
    }

    // Copy edges.
    const auto& new_edges_on_plane = rNewMesh.GetEdgesOnPlanes();
    for( IndexType plane_index = 0; plane_index < 6; ++plane_index){
        const auto& edges = new_edges_on_plane[plane_index];
        for( auto& edge : edges ){
            rTriangleMesh.AddEdgeOnPlane(plane_index, std::get<0>(edge)+initial_number_vertices,
                                                      std::get<1>(edge)+initial_number_vertices,
                                                      std::get<2>(edge)+initial_number_triangles );
        }
    }
}

Unique<TriangleMesh> MeshUtilities::pGetCuboid(const PointType& rLowerPoint, const PointType& rUpperPoint){
    //
    //     2_______3                 y
    //     /      /|                ´|`
    //   6/_____7/ |                 |-->x
    //    | 0   |  /1               /
    //    |     | /                Z
    //   4|____5|/
    //
    auto p_new_triangle_mesh = MakeUnique<TriangleMesh>();

    p_new_triangle_mesh->AddVertex( {rLowerPoint[0], rLowerPoint[1], rLowerPoint[2]} ); //0
    p_new_triangle_mesh->AddVertex( {rUpperPoint[0], rLowerPoint[1], rLowerPoint[2]} ); //1
    p_new_triangle_mesh->AddVertex( {rLowerPoint[0], rUpperPoint[1], rLowerPoint[2]} ); //2
    p_new_triangle_mesh->AddVertex( {rUpperPoint[0], rUpperPoint[1], rLowerPoint[2]} ); //3
    p_new_triangle_mesh->AddVertex( {rLowerPoint[0], rLowerPoint[1], rUpperPoint[2]} ); //4
    p_new_triangle_mesh->AddVertex( {rUpperPoint[0], rLowerPoint[1], rUpperPoint[2]} ); //5
    p_new_triangle_mesh->AddVertex( {rLowerPoint[0], rUpperPoint[1], rUpperPoint[2]} ); //6
    p_new_triangle_mesh->AddVertex( {rUpperPoint[0], rUpperPoint[1], rUpperPoint[2]} ); //7

    // negative x
    p_new_triangle_mesh->AddTriangle({0, 6, 2});
    p_new_triangle_mesh->AddNormal({-1.0, 0.0, 0.0});
    p_new_triangle_mesh->AddTriangle({0, 4, 6});
    p_new_triangle_mesh->AddNormal({-1.0, 0.0, 0.0});

    // postive x
    p_new_triangle_mesh->AddTriangle({1, 7, 5});
    p_new_triangle_mesh->AddNormal({1.0, 0.0, 0.0});
    p_new_triangle_mesh->AddTriangle({1, 3, 7});
    p_new_triangle_mesh->AddNormal({1.0, 0.0, 0.0});

    // negative y
    p_new_triangle_mesh->AddTriangle({4, 1, 5});
    p_new_triangle_mesh->AddNormal({0.0, -1.0, 0.0});
    p_new_triangle_mesh->AddTriangle({4, 0, 1});
    p_new_triangle_mesh->AddNormal({0.0, -1.0, 0.0});

    // postive y
    p_new_triangle_mesh->AddTriangle({6, 7, 3});
    p_new_triangle_mesh->AddNormal({0.0, 1.0, 0.0});
    p_new_triangle_mesh->AddTriangle({6, 3, 2});
    p_new_triangle_mesh->AddNormal({0.0, 1.0, 0.0});

    // negative z
    p_new_triangle_mesh->AddTriangle({1, 0, 3});
    p_new_triangle_mesh->AddNormal({0.0, 0.0, -1.0});
    p_new_triangle_mesh->AddTriangle({0, 2, 3});
    p_new_triangle_mesh->AddNormal({0.0, 0.0, -1.0});

    // positive z
    p_new_triangle_mesh->AddTriangle({4, 5, 7});
    p_new_triangle_mesh->AddNormal({0.0, 0.0, 1.0});
    p_new_triangle_mesh->AddTriangle({4, 7, 6});
    p_new_triangle_mesh->AddNormal({0.0, 0.0, 1.0});

    return std::move(p_new_triangle_mesh);
}

double MeshUtilities::Volume(const TriangleMesh& rTriangleMesh){
    double volume = 0.0;
    const IndexType num_triangles = rTriangleMesh.NumOfTriangles();
    // Loop over all triangles
    for( IndexType i = 0; i < rTriangleMesh.NumOfTriangles(); ++i ){
        const auto p_points = rTriangleMesh.pGetIPsGlobal(i, 3);
        const auto r_points = *p_points;
        // Loop over all points.
        for( const auto& point : r_points ){
            const auto& normal = point.Normal();
            double integrand = Math::Dot(normal, point);
            volume += 1.0/3.0*integrand * point.GetWeight();
        }
    }
    return volume;
}

double MeshUtilities::VolumeOMP(const TriangleMesh& rTriangleMesh){
    double volume = 0.0;
    const IndexType num_triangles = rTriangleMesh.NumOfTriangles();
    // Loop over all triangles in omp parallel.
    #pragma omp parallel for reduction(+ : volume)
    for( IndexType i = 0; i < num_triangles; ++i ){
        const auto p_points = rTriangleMesh.pGetIPsGlobal(i, 3);
        const auto r_points = *p_points;
        // Loop over all points.
        for( const auto& point : r_points ){
            const auto& normal = point.Normal();
            double integrand = Math::Dot(normal, point);
            volume += 1.0/3.0*integrand * point.GetWeight();
        }
    }
    return volume;
}

bool MeshUtilities::IsClosed(const TriangleMesh& rTriangleMesh){
    PointType directional_areas = {0.0, 0.0, 0.0};
    const IndexType num_triangles = rTriangleMesh.NumOfTriangles();
    double total_area = 0.0;
    for( IndexType i = 0; i < num_triangles; ++i ){
        const auto p_points = rTriangleMesh.pGetIPsGlobal(i, 0);
        const auto& r_points = *p_points;
        // Loop over all points.
        for( const auto& point : r_points ){
            const auto& normal = point.Normal();
            directional_areas += normal * point.GetWeight();
            total_area += point.GetWeight();
        }
    }
    double epsilon = std::abs(1.0/total_area*
        (directional_areas[0]+directional_areas[1]+directional_areas[2]));
    return epsilon < 1e-5;
}

} // End namespace tibra