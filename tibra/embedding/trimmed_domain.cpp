// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <random>
///// Project includes
#include "embedding/trimmed_domain.h"
#include "embedding/ray_aabb_primitive.h"
#include "embedding/trimmed_domain_on_plane.h"
#include "utilities/mesh_utilities.h"

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

        if( potential_intersections.size() == 0){
            throw std::runtime_error("TrimmedDomain :: IsInsideTrimmedDomain :: No potential intersections found.");
        }
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
            bool back_facing, parallel;
            if( ray.intersect(p1, p2, p3, t, u, v, back_facing, parallel) ) {
                if( parallel ){
                    is_on_boundary = true;
                    break;
                }
                double sum_u_v = u+v;
                if( t < EPS2 ){ // origin lies on boundary
                    return false;
                }
                // Ray shoots through boundary.
                if( u < 0.0+EPS3 || v < 0.0+EPS3 || sum_u_v > 1.0-EPS3 ){
                    is_on_boundary = true;
                    break;
                }
                // Get orientation of closest triangle. If backfacing=true, triangle is inside.
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
    auto p_trimmed_domain_upper_x = std::make_unique<TrimmedDomainOnPlane>(0, upper_bound, mLowerBound, mUpperBound, this);
    auto p_trimmed_domain_upper_y = std::make_unique<TrimmedDomainOnPlane>(1, upper_bound, mLowerBound, mUpperBound, this);
    auto p_trimmed_domain_upper_z = std::make_unique<TrimmedDomainOnPlane>(2, upper_bound, mLowerBound, mUpperBound, this);
    // Construct trimmed domain on plane z-lower bound of AABB.
    upper_bound = false;
    auto p_trimmed_domain_lower_x = std::make_unique<TrimmedDomainOnPlane>(0, upper_bound, mLowerBound, mUpperBound, this);
    auto p_trimmed_domain_lower_y = std::make_unique<TrimmedDomainOnPlane>(1, upper_bound, mLowerBound, mUpperBound, this);
    auto p_trimmed_domain_lower_z = std::make_unique<TrimmedDomainOnPlane>(2, upper_bound, mLowerBound, mUpperBound, this);

    if( mpTriangleMesh->NumOfTriangles() > 0 ){

        auto p_t1 = p_trimmed_domain_lower_x->pGetTriangulation( *(mpTriangleMesh.get()) );
        auto p_t2 = p_trimmed_domain_upper_x->pGetTriangulation( *(mpTriangleMesh.get()) );
        auto p_t3 = p_trimmed_domain_lower_y->pGetTriangulation( *(mpTriangleMesh.get()) );
        auto p_t4 = p_trimmed_domain_upper_y->pGetTriangulation( *(mpTriangleMesh.get()) );
        auto p_t5 = p_trimmed_domain_lower_z->pGetTriangulation( *(mpTriangleMesh.get()) );
        auto p_t6 = p_trimmed_domain_upper_z->pGetTriangulation( *(mpTriangleMesh.get()) );

        const IndexType num_triangles = p_t1->NumOfTriangles() + p_t2->NumOfTriangles() + p_t3->NumOfTriangles()
            + p_t4->NumOfTriangles() + p_t5->NumOfTriangles() + p_t6->NumOfTriangles();

        mpTriangleMesh->Reserve(2UL*num_triangles);

        MeshUtilities::Append(*(mpTriangleMesh.get()), *(p_t1.get()));
        MeshUtilities::Append(*(mpTriangleMesh.get()), *(p_t2.get()));
        MeshUtilities::Append(*(mpTriangleMesh.get()), *(p_t3.get()));
        MeshUtilities::Append(*(mpTriangleMesh.get()), *(p_t4.get()));
        MeshUtilities::Append(*(mpTriangleMesh.get()), *(p_t5.get()));
        MeshUtilities::Append(*(mpTriangleMesh.get()), *(p_t6.get()));

        MeshUtilities::Refine(*(mpTriangleMesh.get()), mParameters.MinimumNumberOfTriangles());
        p_boundary_ips->reserve(mpTriangleMesh->NumOfTriangles()*6UL);
        for( IndexType triangle_id = 0; triangle_id < mpTriangleMesh->NumOfTriangles(); ++triangle_id ){
            IndexType method = 3; // Creates 6 points per triangle.
            auto p_new_points = mpTriangleMesh->GetIPsGlobal(triangle_id, method);
            p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
        }

        return std::move(p_boundary_ips);
    } else {
        return nullptr;
    }
}

} // End namespace tibra