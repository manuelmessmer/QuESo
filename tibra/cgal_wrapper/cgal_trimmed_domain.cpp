// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// Project includes
#include "cgal_wrapper/cgal_trimmed_domain.h"

namespace cgal {

typedef TrimmedDomainBase::BoundaryIPVectorPtrType BoundaryIPVectorPtrType;
typedef TrimmedDomainBase::BoundingBox BoundingBox;

bool CGALTrimmedDomain::IsInsideTrimmedDomain(const PointType& rPoint) const{
    //  Create CGAL point.
    CGALPointType point(rPoint[0], rPoint[1], rPoint[2]);
    // Perform actual test.
    CGAL::Bounded_side res = (*mpCGALInsideTest)(point);
    if (res == CGAL::ON_BOUNDED_SIDE) { return true; }
    if (res == CGAL::ON_BOUNDARY) { return false; }

    return false;
}

const BoundingBox CGALTrimmedDomain::GetBoundingBoxOfTrimmedDomain() const{
    // Note: std::numeric_limits<double>::min() returns smallest positive number.
    constexpr double min_limit = std::numeric_limits<double>::lowest();
    constexpr double max_limit = std::numeric_limits<double>::max();

    // Initialize bounding box
    BoundingBox bounding_box = { {max_limit, max_limit, max_limit},
                                 {min_limit, min_limit, min_limit} };

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

BoundaryIPVectorPtrType CGALTrimmedDomain::GetBoundaryIps() const {

    /// Create boundary IPs from intersection mesh.
    auto p_boundary_ips = std::make_unique<BoundaryIPVectorType>();
    const SizeType num_triangles = mpTriangleMesh->NumOfTriangles();
    p_boundary_ips->reserve(3*num_triangles);

    for( IndexType triangle_id = 0; triangle_id < num_triangles; ++triangle_id ) {
            IndexType method = 1; // This will create 3 points per triangle.
            auto p_new_points = mpTriangleMesh->GetIPsGlobal(triangle_id, method);
            p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
    }

    return std::move(p_boundary_ips);
}

} // End namespace cgal