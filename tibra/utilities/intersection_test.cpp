// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// CGAl includes
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

// Projecet includes
#include "utilities/intersection_test.h"
#include "utilities/mapping_utilities.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment;
typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
// typedef CGAL::Face_filtered_graph<SurfaceMeshType> Filtered_graph;
// typedef typename boost::property_map<SurfaceMeshType, CGAL::vertex_point_t>::const_type VPMap;
// typedef typename boost::property_map_value<SurfaceMeshType, CGAL::vertex_point_t>::type Point_3;
// typedef typename boost::graph_traits<SurfaceMeshType>::vertex_descriptor   vertex_descriptor;
// //typedef typename boost::graph_traits<SurfaceMeshType>::vertex_descriptor   face_descriptor;
// typedef typename boost::graph_traits<SurfaceMeshType>::halfedge_descriptor halfedge_descriptor;
// typedef boost::property_map<CGAL::Surface_mesh<Point_3>, CGAL::edge_is_feature_t>::type EIFMap;

// namespace PMP = CGAL::Polygon_mesh_processing;

bool IntersectionTest::IsInsideLocalCoordinates(const PointType& point) const{
    const PointType tmp_point = MappingUtilities::FromLocalToGlobalSpace(point, mLowerPoint, mUpperPoint);
    const Point_3 point_3(tmp_point[0], tmp_point[1], tmp_point[2]);

    return IntersectionTest::IsInside(point_3);
}

bool IntersectionTest::IsInside(const PointType& point) const{
    const Point_3 point_3(point[0], point[1], point[2]);

    return IntersectionTest::IsInside(point_3);
}

bool IntersectionTest::IsInside(const Point_3& point) const{
    // TODO: Add more test points, not only vertices of cube.
    CGAL::Bounded_side res = (mpInsideTest)(point);
    if (res == CGAL::ON_BOUNDED_SIDE) { return true; }
    if (res == CGAL::ON_BOUNDARY) { return false; }

    return false;
}

IntersectionTest::IntersectionStatus IntersectionTest::CheckIntersection(const SurfaceMeshType& rSurfaceMesh, const SurfaceMeshType& rCubeMesh, const PointType& rLowerBound, const PointType& rUpperBound) const{
    // Categorize knot spans based on the 8 element vertices
    auto status = CheckInertsectionViaElementVertices( rLowerBound, rUpperBound );

    // CheckInertsectionViaElementVertices does not neccessarily find all intersections
    if( status != IntersectionStatus::Trimmed ){
        const Point_3 point1(rLowerBound[0], rLowerBound[1], rLowerBound[2]);
        const Point_3 point2(rUpperBound[0], rUpperBound[1], rUpperBound[2]);
        const CGAL::Iso_cuboid_3<K> tmp_cuboid( point1, point2, 0);
        if( mAABBTree.do_intersect(tmp_cuboid) ){ // Perform inexact (conservative), but fast test based on AABBTree
            if( CGAL::Polygon_mesh_processing::do_intersect(rSurfaceMesh, rCubeMesh) ){ // Perform exact check
                status = IntersectionStatus::Trimmed;
            }
        }
    }

    return status;
}

IntersectionTest::IntersectionStatus IntersectionTest::CheckIntersection(const SurfaceMeshType& rSurfaceMesh, const SurfaceMeshType& rCubeMesh, const Element& rElement ) const{
    // Categorize knot spans based on the 8 element vertices
    const auto lower_point = rElement.GetGlobalLowerPoint();
    const auto upper_point = rElement.GetGlobalUpperPoint();

    return CheckIntersection(rSurfaceMesh, rCubeMesh, lower_point, upper_point);
}



IntersectionTest::IntersectionStatus IntersectionTest::CheckInertsectionViaElementVertices(const PointType& rLowerBound, const PointType& rUpperBound ) const {

    const PointType lower_point = rLowerBound;
    const PointType upper_point = rUpperBound;

    const PointType point_1 = {upper_point[0], lower_point[1], lower_point[2]};
    const PointType point_2 = {lower_point[0], lower_point[1], upper_point[2]};
    const PointType point_3 = {lower_point[0], lower_point[1], lower_point[2]};
    const PointType point_4 = {lower_point[0], upper_point[1], lower_point[2]};
    const PointType point_5 = {upper_point[0], lower_point[1], upper_point[2]};
    const PointType point_6 = {lower_point[0], upper_point[1], upper_point[2]};
    const PointType point_7 = {upper_point[0], upper_point[1], lower_point[2]};
    const PointType point_8 = {upper_point[0], upper_point[1], upper_point[2]};

    const std::array<PointType,8> points = {point_1, point_2, point_3, point_4, point_5, point_6, point_7, point_8};

    int nb_inside = 0;
    for( auto point : points){
        if( IsInside(point) ) { ++nb_inside; }
    }
    IntersectionTest::IntersectionStatus status;
    if(nb_inside == 0)
        status = Outside;
    else if(nb_inside == 8)
        status = Inside;
    else
        status = Trimmed;

    return status;
}

// bool IntersectionTest::GetIntersectionMesh(const SurfaceMeshType& rSurfaceMesh, SurfaceMeshType& rCubeMesh, Element& rElement, const Parameters& rParam) {

//     // Get all triangles that potentially intersect the element
//     // Conservative and fast search via AABB tree.
//     auto segmented_mesh = GetAllTrianglesThatIntersect(rSurfaceMesh, rCubeMesh, rElement);
//     //std::cout << "000000000000000000000" << std::endl;
//     // Compute the clipped boundary surface that is inside the element.
//     SurfaceMeshType clipped_boundary;
//     CGAL::copy_face_graph(*segmented_mesh, clipped_boundary);
//     // Todo att exception handling
//     //std::cout << "1111111111111111111111" << std::endl;
//     bool success = false;
//     try {
//         //std::cout << "2222222222222222" << std::endl;
//         success = CGAL::Polygon_mesh_processing::clip(clipped_boundary, rCubeMesh, CGAL::Polygon_mesh_processing::parameters::clip_volume(false).throw_on_self_intersection(true));
//         //std::cout << "33333333333333333333" << std::endl;
//     }
//     catch(const std::exception& exc) {
//         std::cout << "Clip not succesful: " << exc.what() << std::endl;
//         return 0;
//     }

//     if(!success){
//         std::cout << "Clip not succesful" << std::endl;
//         return 0;
//     }

//     // Note rCubeMesh is corefinent. Vertices and Edges are added at the intersection with clipped_boundary.
//     // In the next step, we remove all faces that are connected to a corner vertex that is outside the material domain.
//     VPMap vpmap = CGAL::get(CGAL::vertex_point, rCubeMesh);
//     std::vector<vertex_descriptor> vids;
//     const auto point_begin_it = rCubeMesh.vertices_begin();
//     // Loop over first 8 points. These are the 8 corner vertices.
//     //std::cout << "44444444444444444444444" << std::endl;
//     for(IndexType i = 0; i < 8; ++i)
//     {
//         auto point_it = point_begin_it + i;
//         const Point_3& p = get(vpmap, *point_it);
//         if( !IsInside(p) ){ // If outside
//             vids.push_back(*point_it);
//         }
//     }
//     //std::cout << "555555555555" << std::endl;
//     // Mark faces that are connected to the elements in vids.
//     std::list<Primitive_id> selected_faces{};
//     FCCmap fccmap_cube = rCubeMesh.add_property_map<face_descriptor_type, std::size_t>("f:sid").first;
//     for(face_descriptor_type f : faces(rCubeMesh))
//     {
//         bool select = true;
//         for(halfedge_descriptor h : halfedges_around_face(halfedge(f, rCubeMesh), rCubeMesh))
//         {
//             if(std::find(vids.begin(), vids.end(), target(h, rCubeMesh)) != vids.end()){
//                 select = false; }
//         }
//         if( select ){
//             selected_faces.push_back(f);
//         }
//     }
//     //std::cout << "6666666666666" << std::endl;
//     // Get respective mesh, that represent one half of the boundary of the trimmed domain
//     Filtered_graph segment_mesh_graph(rCubeMesh, 0, fccmap_cube);
//     //std::cout << "676767" << std::endl;
//     segment_mesh_graph.set_selected_faces(selected_faces);
//     //std::cout << "6876868" << std::endl;
//     //std::cout << "num faces: " << segment_mesh_graph.number_of_faces() << std::endl;
//     SurfaceMeshType clipped_cube;
//     CGAL::copy_face_graph(segment_mesh_graph, clipped_cube);
//     //std::cout << "6876869" << std::endl;
//     //IO::WriteMeshToVTK(clipped_cube, "cccc.vtk", true);
//     //IO::WriteMeshToVTK(clipped_boundary, "cccc2.vtk", true);
//     //std::cout << "Valid: " << segment_mesh_graph.is_selection_valid() << std::endl;


//     CGAL::copy_face_graph(clipped_boundary, clipped_cube);
//     //std::cout << "77777777777" << std::endl;
//     // Combine both meshes
//     CGAL::PMP::stitch_boundary_cycles(clipped_cube);
//     //CGAL::PMP::remove_self_intersections(clipped_cube);

//     //std::cout << "888888888888" << std::endl;
//     if( !CGAL::PMP::does_bound_a_volume(clipped_cube) ){
//         std::cout << "Mesh Utilities :: Clipped boundary does not bound a volume" << std::endl;
//         //IO::WriteMeshToVTK(clipped_cube, "not_sucess.vtk", false);
//         //return 0;
//     }
//     //std::cout << "9999999999999" << std::endl;
//     if( !CGAL::is_closed(clipped_cube) ){
//         std::cout << "Mesh Utilities :: not closed a volume" << std::endl;
//         // Does almost work, but no always closed mesh
//     }

//     success = Remesh(clipped_cube, rParam);
//     //std::cout << "101010101010" << std::endl;
//     if( !CGAL::PMP::does_bound_a_volume(clipped_cube) ){
//         std::cout << "Mesh Utilities :: not bound a volume" << std::endl;
//     }


//     if( success ){


//         std::cout << "Num faces: " << clipped_cube.num_faces() << std::endl;
//         IO::WriteMeshToVTK(clipped_cube, "remeshed.vtk", false);
//         auto tmp_mesh_ptr = std::make_unique<SurfaceMeshType>(clipped_cube);
//         rElement.pSetSurfaceMesh( tmp_mesh_ptr ); // TODO: Make this move
//         return 1;
//     }


//     std::cout << "no sucess remesh " << std::endl;
//     return 0;
// }


// bool IntersectionTest::Remesh( SurfaceMeshType& rSurfaceMesh, const Parameters& rParam){

//     CGAL::Polygon_mesh_processing::triangulate_faces(rSurfaceMesh);

//     SurfaceMeshType refinend_intersection_mesh{};
//     double edge_length = rParam.InitialTriangleEdgeLength();
//     int iteration_count = 0;

//     while (refinend_intersection_mesh.number_of_faces() < rParam.MinimumNumberOfTriangles() && iteration_count < 10){
//         refinend_intersection_mesh.clear();
//         CGAL::copy_face_graph(rSurfaceMesh, refinend_intersection_mesh);

//         // Element::PositionType positions = refinend_intersection_mesh->points();
//         // for (auto vi = refinend_intersection_mesh->vertices_begin(); vi != refinend_intersection_mesh->vertices_end(); ++vi)
//         // {
//         //     double x = (positions[*vi].x() - rParam.PointA()[0])/std::abs(rParam.PointA()[0] - rParam.PointB()[0]);
//         //     double y = (positions[*vi].y() - rParam.PointA()[1])/std::abs(rParam.PointA()[1] - rParam.PointB()[1]);
//         //     double z = (positions[*vi].z() - rParam.PointA()[2])/std::abs(rParam.PointA()[2] - rParam.PointB()[2]);
//         //     Point_3 p(x,y,z);
//         //     positions[*vi] = p;
//         // }

//         EIFMap eif = CGAL::get(CGAL::edge_is_feature, refinend_intersection_mesh);
//         PMP::detect_sharp_edges(refinend_intersection_mesh, 10, eif);

//         unsigned int nb_iterations = 3;
//         try {
//         CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(refinend_intersection_mesh), edge_length,
//                         refinend_intersection_mesh,  PMP::parameters::number_of_iterations(nb_iterations).use_safety_constraints(false) // authorize all moves
//                                             .edge_is_constrained_map(eif).number_of_relaxation_steps(1));
//         }
//         catch(const std::exception& exc) {
//             if( !CGAL::is_closed(refinend_intersection_mesh) ){
//             std::cout << "Remeshing Exception: " << exc.what() << std::endl;
//             std::cout << "Mesh is not closed. Knot spans will be neglected" << std::endl;
//             return 0;
//             }
//         }
//         edge_length = edge_length * 0.95*std::sqrt((double)refinend_intersection_mesh.number_of_faces() / (double) rParam.MinimumNumberOfTriangles()); // Todo: Make this better!!!
//         iteration_count++;
//     }

//     if (refinend_intersection_mesh.number_of_faces() < rParam.MinimumNumberOfTriangles() ){
//         std::cout << "Warning:: Targeted number of triangles is not reached: "
//         << refinend_intersection_mesh.number_of_faces() << " / " <<  rParam.MinimumNumberOfTriangles() << std::endl;
//     }

//     rSurfaceMesh = refinend_intersection_mesh;
//     // rSurfaceMesh.clear();
//     // CGAL::copy_face_graph(refinend_intersection_mesh, rSurfaceMesh);

//     return 1;
// }

// IntersectionTest::SurfaceMeshPtrType IntersectionTest::GetAllTrianglesThatIntersect(const SurfaceMeshType& rSurfaceMesh, const SurfaceMeshType& rCubeMesh, const Element& rElement) const {

//     const auto lower_point = rElement.GetGlobalLowerPoint();
//     const auto upper_point = rElement.GetGlobalUpperPoint();
//     const Point_3 point1(lower_point[0], lower_point[1], lower_point[2]);
//     const Point_3 point2(upper_point[0], upper_point[1], upper_point[2]);
//     const CGAL::Iso_cuboid_3<K> tmp_cuboid( point1, point2, 0);

//     std::list<Primitive_id> intersections;
//     mAABBTree.all_intersected_primitives(tmp_cuboid, std::back_inserter(intersections));

//     Filtered_graph segment_mesh_graph(rSurfaceMesh, 0, mFCCmap);
//     segment_mesh_graph.set_selected_faces(intersections);

//     SurfaceMeshType tmp_mesh;
//     CGAL::copy_face_graph(segment_mesh_graph, tmp_mesh);

//     return std::make_unique<SurfaceMeshType>(tmp_mesh);
// }