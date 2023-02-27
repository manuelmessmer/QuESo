// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <random>
//// Project includes
#include "embedding/brep_operator.h"
#include "embedding/ray_aabb_primitive.h"

namespace tibra {

typedef BRepOperator::TrimmedDomainBasePtrType TrimmedDomainBasePtrType;

bool BRepOperator::IsInside(const PointType& rPoint) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> drandon(0.5, 1.5);

    if( mTree.IsWithinBoundingBox(rPoint)) {
        bool is_on_boundary = true;
        int intersection_count = 0;
        IndexType iteration = 0UL;
        IndexType max_iteration = 100UL;
        while( is_on_boundary ){
            if( iteration >= max_iteration){ return false; }
            iteration++;
            // Get random direction. Must be postive! -> x>0, y>0, z>0
            Vector3d direction = {drandon(gen), drandon(gen), drandon(gen)};

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
            intersection_count = 0;
            is_on_boundary = false;
            for( auto r : potential_intersections){
                const auto& p1 = mTriangleMesh.P1(r);
                const auto& p2 = mTriangleMesh.P2(r);
                const auto& p3 = mTriangleMesh.P3(r);
                double t, u, v;
                bool back_facing, parallel;
                // Test for actual intersection, if area is area is not zero.
                if( mTriangleMesh.Area(r) > 100.0*ZEROTOL
                        && ray.intersect(p1, p2, p3, t, u, v, back_facing, parallel) ) {
                    intersection_count++;
                    const double sum_u_v = u+v;
                    if( t < ZEROTOL ){ // Origin lies on boundary
                        return false;
                    }
                    // Ray shoots through boundary of triangle.
                    if( u < 0.0+ZEROTOL || v < 0.0+ZEROTOL || sum_u_v > 1.0-ZEROTOL ){
                        is_on_boundary = true;
                        break;
                    }
                    if( parallel ){ // Triangle is parallel to ray.
                        is_on_boundary = true;
                        break;
                    }
                }
            }
        }
        if( intersection_count % 2 == 1){ // If intersection count is odd, return true.
            return true;
        }
    }

    return false;
}

IntersectionStatus BRepOperator::GetIntersectionState(
        const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance) const
{
    // Test if center is inside or outside.
    const PointType center = { 0.5*(rUpperBound[0]+rLowerBound[0]), 0.5*(rUpperBound[1]+rLowerBound[1]), 0.5*(rUpperBound[2]+rLowerBound[2]) };

    IntersectionStatus status = (IsInside(center)) ? Inside : Outside;

    // IntersectionStatus status_confirm = (IsInside(center)) ? Inside : Outside;
    // while( status != status_confirm){
    //     status = (IsInside(center)) ? Inside : Outside;
    //     status_confirm = (IsInside(center)) ? Inside : Outside;
    // }

    AABB_primitive aabb(rLowerBound, rUpperBound);
    auto result = mTree.Query(aabb);

    double snap_tolerance = RelativeSnapTolerance(rLowerBound, rUpperBound, Tolerance);
    for( auto r : result){
        const auto& p1 = mTriangleMesh.P1(r);
        const auto& p2 = mTriangleMesh.P2(r);
        const auto& p3 = mTriangleMesh.P3(r);
        if( aabb.intersect(p1, p2, p3, snap_tolerance) ){
            return IntersectionStatus::Trimmed;
        }
    }

    // If triangle is not intersected, center location will determine if inside or outside.
    return status;
}

TrimmedDomainBasePtrType BRepOperator::GetTrimmedDomain(const PointType& rLowerBound, const PointType& rUpperBound ) const {
    // Instantiate random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> drandon(1, 100);

    // Make copy of lower and upper bound.
    auto lower_bound = rLowerBound;
    auto upper_bound = rUpperBound;

    const double snap_tolerance = RelativeSnapTolerance(lower_bound, upper_bound);
    const double min_volume_ratio = mParameters.Get<double>("min_element_volume_ratio");
    const auto delta = (upper_bound - lower_bound);
    const double volume_non_trimmed_domain = delta[0]*delta[1]*delta[2];

    Unique<TrimmedDomain> best_prev_solution = nullptr;
    double best_error = 1e-2;
    bool switch_plane_orientation = false;
    IndexType iteration = 1UL;
    while( iteration < 10){
        auto p_new_mesh = ClipTriangleMesh(lower_bound, upper_bound);
        if( p_new_mesh->NumOfTriangles() > 0) {
            // Change Direction every intermediate step.
            auto p_trimmed_domain = MakeUnique<TrimmedDomain>(std::move(p_new_mesh), lower_bound, upper_bound, this, mParameters, switch_plane_orientation);
            const auto& r_trimmed_domain_mesh = p_trimmed_domain->GetTriangleMesh();

            double volume = MeshUtilities::Volume( r_trimmed_domain_mesh);
            bool volume_ratio = volume / volume_non_trimmed_domain > min_volume_ratio;

            const double epsilon = MeshUtilities::EstimateQuality(r_trimmed_domain_mesh);
            if( epsilon < 1e-5  ){
                if( volume_ratio )
                    return std::move(p_trimmed_domain);
                else
                    return nullptr;
            }
            else {
                if( epsilon < best_error && volume_ratio ){
                    best_error = epsilon;
                    best_prev_solution = std::move(p_trimmed_domain);
                }
                if( switch_plane_orientation ){
                    PointType lower_perturbation = {drandon(gen)*snap_tolerance, drandon(gen)*snap_tolerance, drandon(gen)*snap_tolerance};
                    PointType upper_perturbation = {drandon(gen)*snap_tolerance, drandon(gen)*snap_tolerance, drandon(gen)*snap_tolerance};

                    lower_bound -= lower_perturbation;
                    upper_bound += upper_perturbation;
                }
            }
        }
        if(switch_plane_orientation ){
            ++iteration;
            switch_plane_orientation = false;
        }
        else {
            switch_plane_orientation = true;
        }
    }

    return best_prev_solution;
}


Unique<std::vector<IndexType>> BRepOperator::GetIntersectedTriangleIds(
        const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance ) const {

    // Perform fast search based on aabb tree. Conservative search.
    AABB_primitive aabb(rLowerBound, rUpperBound);
    auto potential_intersections = mTree.Query(aabb);

    auto intersected_triangle_ids = MakeUnique<std::vector<IndexType>>();
    int count = 0;
    intersected_triangle_ids->reserve(potential_intersections.size());
    for( auto triangle_id : potential_intersections){
        const auto& p1 = mTriangleMesh.P1(triangle_id);
        const auto& p2 = mTriangleMesh.P2(triangle_id);
        const auto& p3 = mTriangleMesh.P3(triangle_id);

        // Perform actual intersection test.
        if( aabb.intersect(p1, p2, p3, Tolerance) ){
            intersected_triangle_ids->push_back(triangle_id);
        }
    }

    return intersected_triangle_ids;
}

Unique<TriangleMesh> BRepOperator::ClipTriangleMesh(
        const PointType& rLowerBound, const PointType& rUpperBound ) const {

    const double snap_tolerance = 0.1*RelativeSnapTolerance(rUpperBound, rLowerBound);
    auto p_intersected_triangle_ids = GetIntersectedTriangleIds(rLowerBound, rUpperBound, snap_tolerance);
    auto p_triangle_mesh = MakeUnique<TriangleMesh>();
    p_triangle_mesh->Reserve( 2*p_intersected_triangle_ids->size() );
    p_triangle_mesh->ReserveEdgesOnPlane( p_intersected_triangle_ids->size() );
    for( auto triangle_id : (*p_intersected_triangle_ids) ){
        const auto& P1 = mTriangleMesh.P1(triangle_id);
        const auto& P2 = mTriangleMesh.P2(triangle_id);
        const auto& P3 = mTriangleMesh.P3(triangle_id);
        const auto& normal = mTriangleMesh.Normal(triangle_id);
        auto p_polygon = Clipper::ClipTriangle(P1, P2, P3, normal, rLowerBound, rUpperBound );
        if( p_polygon ){
            p_polygon->AddToTriangleMesh(*p_triangle_mesh.get());
        }
    }
    if( p_triangle_mesh->NumOfTriangles() > 0){
        p_triangle_mesh->Check();
    }
    return std::move(p_triangle_mesh);
}


Unique<TriangleMesh> BRepOperator::ClipTriangleMeshUnique(const PointType& rLowerBound, const PointType& rUpperBound ) const {
    // This needs improvement. Global triangle clipping.
    const PointType offset(30*ZEROTOL, 30*ZEROTOL, 30*ZEROTOL);
    const auto lower_bound = rLowerBound + offset;
    const auto upper_bound = rUpperBound + offset;
    double snap_tolerance = 1.0*ZEROTOL;

    auto p_intersected_triangle_ids = GetIntersectedTriangleIds(lower_bound, upper_bound, snap_tolerance);
    auto p_triangle_mesh = MakeUnique<TriangleMesh>();
    p_triangle_mesh->Reserve( 2*p_intersected_triangle_ids->size() );
    p_triangle_mesh->ReserveEdgesOnPlane( p_intersected_triangle_ids->size() );

    for( auto triangle_id : (*p_intersected_triangle_ids) ){
        const auto& P1 = mTriangleMesh.P1(triangle_id);
        const auto& P2 = mTriangleMesh.P2(triangle_id);
        const auto& P3 = mTriangleMesh.P3(triangle_id);
        const auto& normal = mTriangleMesh.Normal(triangle_id);
        auto p_polygon = Clipper::ClipTriangle(P1, P2, P3, normal, lower_bound, upper_bound);
        if( p_polygon ){
            p_polygon->AddToTriangleMesh(*p_triangle_mesh.get());
        }
    }
    snap_tolerance = RelativeSnapTolerance(rLowerBound, rUpperBound);
    if( p_triangle_mesh->NumOfTriangles() > 0){
        p_triangle_mesh->Check();
    }
    return std::move(p_triangle_mesh);
}

} // End namespace tibra

// Winding numbers algorithm: It actually works!!!

// TIBRA_INFO << "Done" << std::endl;
// TIBRA_INFO << "Tree timer: " << tree_timer << std::endl;
// TIBRA_INFO << "ray timer: " << ray_timer << std::endl;
// return MakeUnique<std::vector<bool>>(rr);
// for( int triangle_id = 0; triangle_id < rTriangleMesh.NumOfTriangles(); ++triangle_id ){
//     const auto& p1 = rTriangleMesh.P1(triangle_id);
//     const auto& p2 = rTriangleMesh.P2(triangle_id);
//     const auto& p3 = rTriangleMesh.P3(triangle_id);

//     IndexType count = 0;
//     for( const auto& point : rPoints){
//     //const auto point = rPoints[0];
//         PointType a{p1[0] - point[0], p1[1] - point[1], p1[2] - point[2]};
//         PointType b{p2[0] - point[0], p2[1] - point[1], p2[2] - point[2]};
//         PointType c{p3[0] - point[0], p3[1] - point[1], p3[2] - point[2]};

//         // Compute determinant of:
//         // [a0 b0 c0]
//         // [a1 b1 c1]
//         // [a2 b2 c2]
//         const long double omega =  a[0] * (b[1]*c[2] - b[2]*c[1])
//                                   -b[0] * (a[1]*c[2] - a[2]*c[1])
//                                   +c[0] * (a[1]*b[2] - a[2]*b[1]);


//         const long double anorm = std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
//         const long double bnorm = std::sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
//         const long double cnorm = std::sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

//         long double k = anorm*bnorm*cnorm;

//         k += cnorm * ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
//         k += anorm * ( b[0]*c[0] + b[1]*c[1] + b[2]*c[2]);
//         k += bnorm * ( a[0]*c[0] + a[1]*c[1] + a[2]*c[2]);

//         ret[count] += std::atan2(omega,k);
//         count++;
//     }

// }

// Unique<std::vector<bool>> p_result = MakeUnique<std::vector<bool>>(rPoints.size(), false);
// auto& r_result = *p_result;
// for( int i = 0; i < ret.size(); ++i ){
//     if( ret[i] >= (2.0*My_PI-EPS2) ) // Due to limited precision of double small treshhold (EPS2) is required.
//         r_result[i] = true;
// }

// return std::move(p_result);

// def is_inside_turbo(triangles, X):
// 	# Compute euclidean norm along axis 1
// 	def anorm2(X):
// 		return numpy.sqrt(numpy.sum(X ** 2, axis = 1))



// 	# Compute 3x3 determinant along axis 1
// 	def adet(X, Y, Z):
// 		ret  = numpy.multiply(numpy.multiply(X[:,0], Y[:,1]), Z[:,2])
// 		ret += numpy.multiply(numpy.multiply(Y[:,0], Z[:,1]), X[:,2])
// 		ret += numpy.multiply(numpy.multiply(Z[:,0], X[:,1]), Y[:,2])
// 		ret -= numpy.multiply(numpy.multiply(Z[:,0], Y[:,1]), X[:,2])
// 		ret -= numpy.multiply(numpy.multiply(Y[:,0], X[:,1]), Z[:,2])
// 		ret -= numpy.multiply(numpy.multiply(X[:,0], Z[:,1]), Y[:,2])
// 		return ret



// 	# One generalized winding number per input vertex
// 	ret = numpy.zeros(X.shape[0], dtype = X.dtype)

// 	# Accumulate generalized winding number for each triangle
// 	for U, V, W in triangles:
// 		A, B, C = U - X, V - X, W - X
// 		omega = adet(A, B, C)

// 		a, b, c = anorm2(A), anorm2(B), anorm2(C)
// 		k  = a * b * c
// 		k += c * numpy.sum(numpy.multiply(A, B), axis = 1)
// 		k += a * numpy.sum(numpy.multiply(B, C), axis = 1)
// 		k += b * numpy.sum(numpy.multiply(C, A), axis = 1)

// 		ret += numpy.arctan2(omega, k)

// 	# Job done
// 	return ret >= 2 * numpy.pi