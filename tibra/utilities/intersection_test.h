// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef INTERSECTION_TEST_H
#define INTERSECTION_TEST_H

// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
// #include <CGAL/Side_of_triangle_mesh.h>
// #include <CGAL/Surface_mesh.h>
// #include <CGAL/boost/graph/Face_filtered_graph.h>
// #include <boost/property_map/property_map.hpp>
// #include <CGAL/Polygon_mesh_processing/connected_components.h>
// #include <CGAL/Polygon_mesh_processing/clip.h>

// Projec includes
#include "geometries/element.h"
// #include "io/io_utilities.h"

// External includes
#include <iostream>
#include <array>
#include <string>
#include <map>

//@brief This class is deprecated and will be removed soon! All related functions are move to BRepOperator.
class IntersectionTest {

public:
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3;
    typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
    typedef std::unique_ptr<SurfaceMeshType> SurfaceMeshPtrType;
    typedef CGAL::Side_of_triangle_mesh<SurfaceMeshType, K> InsideTestType;
    typedef std::array<double,3> PointType;

    typedef CGAL::AABB_face_graph_triangle_primitive<SurfaceMeshType>           AABBFaceGraphsPrimitivesType;
    typedef CGAL::AABB_traits<K, AABBFaceGraphsPrimitivesType>                  AABBFaceGraphsTraits;
    typedef CGAL::AABB_tree<AABBFaceGraphsTraits>                               AABBTreeType;

    // typedef boost::graph_traits<SurfaceMeshType>::face_descriptor               face_descriptor_type;
    // typedef SurfaceMeshType::Property_map<face_descriptor_type, std::size_t>    FCCmap;

    // typedef AABBTreeType::Primitive_id Primitive_id;

    enum IntersectionStatus {Inside, Outside, Trimmed};

    IntersectionTest(const SurfaceMeshType& Polyhedron, const PointType& LowerPoint, const PointType& UpperPoint) :
        mLowerPoint(LowerPoint), mUpperPoint(UpperPoint), mpInsideTest(Polyhedron)
    {
        // Build AABBTree
        CGAL::Polygon_mesh_processing::build_AABB_tree(Polyhedron, mAABBTree);

        //mFCCmap = const_cast<SurfaceMeshType&>(Polyhedron).add_property_map<face_descriptor_type, std::size_t>("f:sid").first;
    }

    bool IsInsideLocalCoordinates(const PointType& point) const;

    bool IsInside(const PointType& point) const;

    bool IsInside(const Point_3& point) const;

    IntersectionStatus CheckIntersection(const SurfaceMeshType& rSurfaceMesh, const SurfaceMeshType& rCubeMesh, const Element& rElement) const;

    IntersectionStatus CheckIntersection(const SurfaceMeshType& rSurfaceMesh, const SurfaceMeshType& rCubeMesh, const PointType& rLowerBound, const PointType& rUpperBound) const;
    // bool GetIntersectionMesh(const SurfaceMeshType& rSurfaceMesh, SurfaceMeshType& rCubeMesh, Element& rElemen, const Parameters& rParam);

    // bool Remesh( SurfaceMeshType& rSurfaceMesh, const Parameters& rParam);

private:
    IntersectionStatus CheckInertsectionViaElementVertices( const PointType& rLowerBound, const PointType& rUpperBound) const;

    // SurfaceMeshPtrType GetAllTrianglesThatIntersect(const SurfaceMeshType& rSurfaceMesh, const SurfaceMeshType& rCubeMesh, const Element& rElement) const;



    InsideTestType mpInsideTest;
    AABBTreeType mAABBTree;
    PointType mLowerPoint;
    PointType mUpperPoint;
    //FCCmap mFCCmap;

};

#endif // INSIDE_TEST_H