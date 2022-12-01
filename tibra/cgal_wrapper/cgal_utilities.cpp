// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// CGAL includes
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
//// Project includes
#include "cgal_wrapper/cgal_utilities.h"

namespace tibra {
namespace cgal {

bool CGALUtilities::CopyMesh(const CGALMeshType& rInputMesh, TriangleMesh& rOutputMesh ){

    typedef typename boost::graph_traits<CGALMeshType>::vertex_descriptor   vertex_descriptor;
    typedef typename boost::property_map<CGALMeshType, CGAL::vertex_point_t>::const_type VPMap;

    VPMap vpmap = CGAL::get(CGAL::vertex_point, rInputMesh);

    const SizeType num_elements = rInputMesh.number_of_faces();
    const SizeType num_points = rInputMesh.number_of_vertices();

    rOutputMesh.Clear();
    rOutputMesh.Reserve(num_elements);

    std::map<vertex_descriptor, IndexType> vids;
    IndexType inum = 0;
    auto raw_points = rInputMesh.points().data()->cartesian_begin();
    for(vertex_descriptor v : vertices(rInputMesh))
    {
        const CGALPointType& p = get(vpmap, v);
        PointType tmp_point = {CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()) };
        rOutputMesh.AddVertex( tmp_point );
        vids[v] = inum++;
    }

    for(auto f : faces(rInputMesh))
    {
        Vector3i vertex_ids;
        IndexType count = 0;
        for(auto h : halfedges_around_face(halfedge(f, rInputMesh), rInputMesh))
        {
            vertex_ids[count] = vids[target(h, rInputMesh)];
            count++;
        }

        auto cgal_normal = CGAL::Polygon_mesh_processing::compute_face_normal(f, rInputMesh);
        PointType normal = {cgal_normal[0], cgal_normal[1], cgal_normal[2]};
        rOutputMesh.AddTriangle(vertex_ids);
        rOutputMesh.AddNormal(normal);
    }

    return 1;
}

bool CGALUtilities::CopyMesh(const TriangleMesh& rInputMesh, CGALMeshType& rOutputMesh ){
    std::map<IndexType, CGALMeshType::Vertex_index> index_map{};
    const auto& r_vertices = rInputMesh.GetVertices();
    const auto v_it_begin = r_vertices.begin();
    for( IndexType i = 0; i < rInputMesh.NumOfVertices(); ++i){
        auto v = *(v_it_begin + i);
        auto index1 = rOutputMesh.add_vertex(CGALPointType(v[0], v[1], v[2]));
        index_map.insert( std::pair<IndexType, CGALMeshType::Vertex_index>( i, index1) );
    }

    for( IndexType i = 0; i < rInputMesh.NumOfTriangles(); ++i){
        const auto& ids = rInputMesh.VertexIds(i);
        rOutputMesh.add_face( index_map[ids[0]], index_map[ids[1]], index_map[ids[2]]  );
    }

    return 1;
}

} // End namespace cgal
} // End namespace tibra