// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <random>
///// Project includes
#include "embedding/trimmed_domain.h"
#include "embedding/ray_aabb_primitive.h"
#include "embedding/trimmed_domain_on_plane.h"

namespace tibra {

typedef TrimmedDomainBase::BoundaryIPVectorPtrType BoundaryIPVectorPtrType;
typedef TrimmedDomainBase::BoundingBox BoundingBox;

bool TrimmedDomain::IsInsideTrimmedDomain(const PointType& rPoint) const {

    const IndexType num_triangles = mpTriangleMesh->NumOfTriangles();
    if( num_triangles == 0){
        return true;
    }
    bool is_on_boundary = true;
    bool is_inside = false;
    int intersection_count = 0;
    IndexType current_id = 0;
    while( is_on_boundary ){
        if( current_id >= num_triangles ){
            return false;
        }
        // Get direction
        const auto center_triangle = mpTriangleMesh->Center(current_id);
        Vector3d direction = center_triangle - rPoint;

        // Normalize
        double norm_direction = direction.Norm();
        direction /= norm_direction;

        // Construct ray
        Ray_AABB_primitive ray(rPoint, direction);
        // Get potential ray intersections from AABB tree.
        auto potential_intersections = mTree.Query(ray);

        // Test if potential intersections actually intersect.
        // If intersection lies on boundary cast a new ray.
        // @todo Use symbolic perturbations: http://dl.acm.org/citation.cfm?id=77639
        is_on_boundary = false;
        double min_distance = std::numeric_limits<double>::max();
        is_inside = false;
        for( auto r : potential_intersections){
            const auto& p1 = mpTriangleMesh->P1(r);
            const auto& p2 = mpTriangleMesh->P2(r);
            const auto& p3 = mpTriangleMesh->P3(r);
            double t, u, v;
            bool back_facing;
            if( ray.intersect(p1, p2, p3, t, u, v, back_facing) ) {
                double sum_u_v = u+v;
                if( t < 1e-10 ){ // origin lies on boundary
                    return false;
                }
                if( u < 0.0+1e-14 || v < 0.0+1e-14 || sum_u_v > 1.0-1e-14 ){
                    is_on_boundary = true; // Hits boundary
                    break;
                }
                if( t < min_distance ){
                    is_inside = back_facing;
                    min_distance = t;
                }
            }
        }
        current_id++;
    }
    return is_inside;
}

const BoundingBox TrimmedDomain::GetBoundingBoxOfTrimmedDomain() const {
    // Note: std::numeric_limits<double>::min() returns smallest positive number.
    constexpr double min_limit = std::numeric_limits<double>::lowest();
    constexpr double max_limit = std::numeric_limits<double>::max();

    // Initialize bounding box
    BoundingBox bounding_box = { {max_limit, max_limit, max_limit},
                                 {min_limit, min_limit, min_limit} };

    // Since TrimmedDomain holds only a clipped TriangleMesh (mesh is not closed), vertices of AABB must also be considered.
    PointType point_1(mUpperBound[0], mLowerBound[1], mLowerBound[2]);
    PointType point_2(mLowerBound[0], mLowerBound[1], mUpperBound[2]);
    PointType point_3(mLowerBound[0], mLowerBound[1], mLowerBound[2]);
    PointType point_4(mLowerBound[0], mUpperBound[1], mLowerBound[2]);
    PointType point_5(mUpperBound[0], mLowerBound[1], mUpperBound[2]);
    PointType point_6(mLowerBound[0], mUpperBound[1], mUpperBound[2]);
    PointType point_7(mUpperBound[0], mUpperBound[1], mLowerBound[2]);
    PointType point_8(mUpperBound[0], mUpperBound[1], mUpperBound[2]);
    std::array<PointType,8> points = {point_1, point_2, point_3, point_4, point_5, point_6, point_7, point_8};

    for( auto& p : points){
        if( IsInsideTrimmedDomain(p) ){
            // Loop over all 3 dimensions
            for( IndexType i = 0; i < 3; ++i){
                if( p[i] < bounding_box.first[i] ){ // Find min values
                    bounding_box.first[i] = p[i];
                }
                if( p[i] > bounding_box.second[i] ){ // Find max values
                    bounding_box.second[i] = p[i];
                }
            }
        }
    }

    const auto& vertices = mpTriangleMesh->GetVertices();
    for( auto& v : vertices ){
        // Loop over all 3 dimensions
        for( IndexType i = 0; i < 3; ++i){
            if( v[i] < bounding_box.first[i] ){ // Find min values
                bounding_box.first[i] = v[i];
            }
            if( v[i] > bounding_box.second[i] ){ // Find max values
                bounding_box.second[i] = v[i];
            }
        }
    }

    return bounding_box;
}

BoundaryIPVectorPtrType TrimmedDomain::pGetBoundaryIps() const{
    // Pointer to boundary integration points
    auto p_boundary_ips = std::make_unique<BoundaryIPVectorType>();

    // Construct trimmed domain on plane z-upper bound of AABB.
    bool upper_bound = true;
    IndexType plane_index = 2;
    auto p_trimmed_domain_upper = std::make_unique<TrimmedDomainOnPlane>(plane_index, upper_bound, mLowerBound, mUpperBound, this);
    // Construct trimmed domain on plane z-lower bound of AABB.
    upper_bound = false;
    auto p_trimmed_domain_lower = std::make_unique<TrimmedDomainOnPlane>(plane_index, upper_bound, mLowerBound, mUpperBound, this);

    mpTriangleMesh->Refine(0.005);
    // Loop over all intersected triangles
    for( IndexType triangle_id = 0; triangle_id < mpTriangleMesh->NumOfTriangles(); ++triangle_id ){
        const auto& P1 = mpTriangleMesh->P1(triangle_id);
        const auto& P2 = mpTriangleMesh->P2(triangle_id);
        const auto& P3 = mpTriangleMesh->P3(triangle_id);
        const auto& normal = mpTriangleMesh->Normal(triangle_id);

        // Triangle is fully contained, does not need to be clipped.
        // @TODO: Add to mesh and enable the option to refine the mesh, if number of triangles are belowe certain treshhold.
        auto p_new_points = mpTriangleMesh->GetIPsGlobal(triangle_id, 1);
        p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());

        /// Upper plane
        if( IsPointOnUpperPlane(P1, plane_index) && IsPointOnUpperPlane(P2, plane_index) ){
            std::array<PointType,2> new_edge = {P1, P2};
            p_trimmed_domain_upper->InsertEdge(new_edge, normal);
        }
        if( IsPointOnUpperPlane(P2, plane_index) && IsPointOnUpperPlane(P3, plane_index) ){
            std::array<PointType,2> new_edge = {P2, P3};
            p_trimmed_domain_upper->InsertEdge(new_edge, normal);
        }
        if( IsPointOnUpperPlane(P3, plane_index) && IsPointOnUpperPlane(P1, plane_index) ){
            std::array<PointType,2> new_edge = {P3, P1};
            p_trimmed_domain_upper->InsertEdge(new_edge, normal);
        }

        /// Lower Plane
        if( IsPointOnLowerPlane(P1, plane_index) && IsPointOnLowerPlane(P2, plane_index) ){
            std::array<PointType,2> new_edge = {P1, P2};
            p_trimmed_domain_lower->InsertEdge(new_edge, normal);
        }
        if( IsPointOnLowerPlane(P2, plane_index) && IsPointOnLowerPlane(P3, plane_index) ){
            std::array<PointType,2> new_edge = {P2, P3};
            p_trimmed_domain_lower->InsertEdge(new_edge, normal);
        }
        if( IsPointOnLowerPlane(P3, plane_index) && IsPointOnLowerPlane(P1, plane_index) ){
            std::array<PointType,2> new_edge = {P3, P1};
            p_trimmed_domain_lower->InsertEdge(new_edge, normal);
        }
    }
    /// Mapping:
    //
    //     a_______b                 y
    //     /      /|                ´|`
    //   c/_____d/ |<-- lower plane  |-->x
    //    |     |  /                /
    //    |     |</-- upper plane  Z
    //    |_____|/
    //
    // pGetBoundaryIPs needs to know if a,b,c,d are inside or outside domain.
    PointType point_a = { mLowerBound[0], mUpperBound[1], mLowerBound[2] };
    bool state_a = IsInsideTrimmedDomain(point_a);
    PointType point_b = { mUpperBound[0], mUpperBound[1], mLowerBound[2] };
    bool state_b = IsInsideTrimmedDomain(point_b);
    PointType point_c = { mLowerBound[0], mUpperBound[1], mUpperBound[2] };
    bool state_c = IsInsideTrimmedDomain(point_c);
    PointType point_d = { mUpperBound[0], mUpperBound[1], mUpperBound[2] };
    bool state_d = IsInsideTrimmedDomain(point_d);

    //@TODO: Implement this as well.
    double egde_cube_lenght_ratio = 1.0/10.0;

    p_trimmed_domain_lower->pGetBoundaryIPs( p_boundary_ips, state_a, state_b );
    p_trimmed_domain_upper->pGetBoundaryIPs( p_boundary_ips, state_c, state_d );

    return std::move(p_boundary_ips);
}

} // End namespace tibra