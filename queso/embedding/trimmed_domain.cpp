// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <random>
//// Project includes
#include "includes/define.hpp"
#include "embedding/trimmed_domain.h"
#include "embedding/ray_aabb_primitive.h"
#include "utilities/mesh_utilities.h"
#include "io/io_utilities.h"

namespace queso {

typedef TrimmedDomainBase::BoundaryIPVectorPtrType BoundaryIPVectorPtrType;
typedef TrimmedDomainBase::BoundingBox BoundingBox;

bool TrimmedDomain::IsInsideTrimmedDomain(const PointType& rPoint, bool& rSuccess) const {

    const IndexType num_triangles = mClippedMesh.NumOfTriangles();
    if( num_triangles == 0){
        return true;
    }
    rSuccess = true;
    bool success_local = false;
    bool is_inside = false;
    IndexType current_id = 0;
    while( !success_local ){
        // Return false if all triangles are tested, but non valid (all are parallel or on_boundary)
        if( current_id >= num_triangles ){
            rSuccess = false;
            return false;
        }

        // Get direction
        const auto center_triangle = mClippedMesh.Center(current_id);
        Vector3d direction = center_triangle - rPoint;

        // Normalize
        double norm_direction = direction.Norm();
        direction /= norm_direction;

        // Construct ray
        Ray_AABB_primitive ray(rPoint, direction);

        // Get vertices of current triangle
        const auto& p1 = mClippedMesh.P1(current_id);
        const auto& p2 = mClippedMesh.P2(current_id);
        const auto& p3 = mClippedMesh.P3(current_id);

        // Make sure target triangle is not parallel and has a significant area.
        const double area = mClippedMesh.Area(current_id);
        if( !ray.is_parallel(p1, p2, p3, 100.0*mSnapTolerance) && area >  100*ZEROTOL) {
            std::tie(is_inside, success_local) = mGeometryQuery.IsInside(ray);
        }
        current_id++;
    }
    return is_inside;
}

const BoundingBox TrimmedDomain::GetBoundingBoxOfTrimmedDomain() const {
    // Initialize bounding box
    BoundingBox bounding_box = { {MAXD, MAXD, MAXD},
                                 {LOWESTD, LOWESTD, LOWESTD} };

    // Loop over all vertices
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
    auto p_boundary_ips = MakeUnique<BoundaryIPVectorType>();

    p_boundary_ips->reserve(mpTriangleMesh->NumOfTriangles()*6UL);
    for( IndexType triangle_id = 0; triangle_id < mpTriangleMesh->NumOfTriangles(); ++triangle_id ){
        IndexType method = 3; // Creates 6 points per triangle.
        auto p_new_points = mpTriangleMesh->pGetIPsGlobal(triangle_id, method);
        p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
    }

    return p_boundary_ips;
}

IntersectionStatusType TrimmedDomain::GetIntersectionState(
        const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance) const
{
    if( mGeometryQuery.DoIntersect(rLowerBound, rUpperBound, Tolerance) ){
        return IntersectionStatus::Trimmed;
    }

    // Test if center is inside or outside.
    const PointType center = (rLowerBound + rUpperBound) * 0.5;
    const auto status = (TrimmedDomainBase::IsInsideTrimmedDomain(center)) ? IntersectionStatus::Inside : IntersectionStatus::Outside;

    // If triangle is not intersected, center location will determine if inside or outside.
    return status;
}

} // End namespace queso