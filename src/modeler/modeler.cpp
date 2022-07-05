// Project includes
#include "modeler/modeler.h"

std::unique_ptr<Modeler::SurfaceMeshType> Modeler::make_cube_3( std::array<double,3> lower_point, std::array<double,3> upper_point) {

    typedef typename Mesh::Point_3         Point;
    typedef typename Mesh::Halfedge_handle Halfedge_handle;

    Mesh P;
    CGAL_precondition( P.is_valid());

    Point point_1(upper_point[0], lower_point[1], lower_point[2]);
    Point point_2(lower_point[0], lower_point[1], upper_point[2]);
    Point point_3(lower_point[0], lower_point[1], lower_point[2]);
    Point point_4(lower_point[0], upper_point[1], lower_point[2]);
    Point point_5(upper_point[0], lower_point[1], upper_point[2]);
    Point point_6(lower_point[0], upper_point[1], upper_point[2]);
    Point point_7(upper_point[0], upper_point[1], lower_point[2]);
    Point point_8(upper_point[0], upper_point[1], upper_point[2]);
    std::array<Point,8> points = {point_1, point_2, point_3, point_4, point_5, point_6, point_7, point_8};

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

    SurfaceMeshType surface_mesh;
    CGAL::copy_face_graph(P, surface_mesh);
    return std::make_unique<SurfaceMeshType>(surface_mesh);
}

void Modeler::make_cube_3( Modeler::SurfaceMeshType& rSurfaceMesh, std::array<double,3> lower_point, std::array<double,3> upper_point) {
    typedef typename Mesh::Point_3           Point;
    typedef typename Mesh::Halfedge_handle   Halfedge_handle;

    Mesh P;
    CGAL_precondition( P.is_valid());

    Point point_1(upper_point[0], lower_point[1], lower_point[2]);
    Point point_2(lower_point[0], lower_point[1], upper_point[2]);
    Point point_3(lower_point[0], lower_point[1], lower_point[2]);
    Point point_4(lower_point[0], upper_point[1], lower_point[2]);
    Point point_5(upper_point[0], lower_point[1], upper_point[2]);
    Point point_6(lower_point[0], upper_point[1], upper_point[2]);
    Point point_7(upper_point[0], upper_point[1], lower_point[2]);
    Point point_8(upper_point[0], upper_point[1], upper_point[2]);
    std::array<Point,8> points = {point_1, point_2, point_3, point_4, point_5, point_6, point_7, point_8};

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
    CGAL::copy_face_graph(P, rSurfaceMesh);
}

vtkSmartPointer<vtkHexahedron> Modeler::GetVTKHexahedron( std::array<double,3> lower_point, std::array<double,3> upper_point){

    auto hexahedron = vtkSmartPointer<vtkHexahedron>::New();
    for (int i = 0; i < hexahedron->GetNumberOfPoints(); ++i)
    {
        hexahedron->GetPointIds()->SetId(i, i);
    }
    hexahedron->GetPoints()->SetPoint(0,lower_point[0], lower_point[1], lower_point[2]);
    hexahedron->GetPoints()->SetPoint(1,upper_point[0], lower_point[1], lower_point[2]);
    hexahedron->GetPoints()->SetPoint(2,upper_point[0], upper_point[1], lower_point[2]);
    hexahedron->GetPoints()->SetPoint(3,lower_point[0], upper_point[1], lower_point[2]);
    hexahedron->GetPoints()->SetPoint(4,lower_point[0], lower_point[1], upper_point[2]);
    hexahedron->GetPoints()->SetPoint(5,upper_point[0], lower_point[1], upper_point[2]);
    hexahedron->GetPoints()->SetPoint(6,upper_point[0], upper_point[1], upper_point[2]);
    hexahedron->GetPoints()->SetPoint(7,lower_point[0], upper_point[1], upper_point[2]);

    return hexahedron;
}