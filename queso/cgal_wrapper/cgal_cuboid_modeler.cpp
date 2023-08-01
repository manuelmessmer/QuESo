// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// Project includes
#include "cgal_wrapper/cgal_cuboid_modeler.h"

namespace queso {
namespace cgal {

typedef CuboidModeler::CGALMeshType CGALMeshType;
typedef CuboidModeler::CGALMeshPtrType CGALMeshPtrType;

CGALMeshPtrType CuboidModeler::MakeCuboid( const PointType& rLowerPoint, const PointType& rUpperPoint) {

    typedef typename CGALPolyhedronMeshType::Halfedge_handle Halfedge_handle;

    CGALPolyhedronMeshType P;
    CGAL_precondition( P.is_valid());

    CGALPointType point_1(rUpperPoint[0], rLowerPoint[1], rLowerPoint[2]);
    CGALPointType point_2(rLowerPoint[0], rLowerPoint[1], rUpperPoint[2]);
    CGALPointType point_3(rLowerPoint[0], rLowerPoint[1], rLowerPoint[2]);
    CGALPointType point_4(rLowerPoint[0], rUpperPoint[1], rLowerPoint[2]);
    CGALPointType point_5(rUpperPoint[0], rLowerPoint[1], rUpperPoint[2]);
    CGALPointType point_6(rLowerPoint[0], rUpperPoint[1], rUpperPoint[2]);
    CGALPointType point_7(rUpperPoint[0], rUpperPoint[1], rLowerPoint[2]);
    CGALPointType point_8(rUpperPoint[0], rUpperPoint[1], rUpperPoint[2]);

    Halfedge_handle h = P.make_tetrahedron( point_1, point_2, point_3, point_4);
    Halfedge_handle g = h->next()->opposite()->next();
    P.split_edge( h->next());
    P.split_edge( g->next());
    P.split_edge( g);
    h->next()->vertex()->point()     = point_5;
    g->next()->vertex()->point()     = point_6;
    g->opposite()->vertex()->point() = point_7;
    Halfedge_handle f = P.split_facet( g->next(),
                                       g->next()->next()->next());
    Halfedge_handle e = P.split_edge( f);
    e->vertex()->point() = point_8;
    P.split_facet( e, f->next()->next());

    CGAL::Polygon_mesh_processing::triangulate_faces(P);
    CGAL_postcondition( P.is_valid());

    CGALMeshType surface_mesh;
    CGAL::copy_face_graph(P, surface_mesh);
    return MakeUnique<CGALMeshType>(surface_mesh);
}

} // End namespace cgal
} // End namespace queso