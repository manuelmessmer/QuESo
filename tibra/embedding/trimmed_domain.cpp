// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <random>
///// Project includes
#include "embedding/trimmed_domain.h"
#include "embedding/ray_aabb_primitive.h"

namespace tibra {



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
        const auto center_triangle = mpTriangleMesh->Center(current_id);
        Vector3d direction = center_triangle - rPoint;

        // Normalize
        double sum_direction = std::sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
        direction /= sum_direction;
        //std::for_each(direction.begin(), direction.end(), [sum_direction](auto& rValue) { rValue /= sum_direction;} );

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

} // End namespace tibra