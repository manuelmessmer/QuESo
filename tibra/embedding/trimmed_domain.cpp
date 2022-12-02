// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <random>
///// Project includes
#include "embedding/trimmed_domain.h"
#include "embedding/ray_aabb_primitive.h"

namespace tibra {

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
        double sum_direction = std::sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
        direction /= sum_direction;

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
                    is_on_boundary = true;
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

    // Check vertices of aabb that are inside trimmed domain;
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

} // End namespace tibra