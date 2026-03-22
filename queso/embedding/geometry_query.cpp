//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

//// Project includes
#include "queso/embedding/geometry_query.h"
#include "queso/utilities/triangle_utilities.hpp"

//// STL includes
#include <optional>

namespace queso {

    bool GeometryQuery::IsWithinBoundingBox(PointView rPoint) const {
        return mTree.IsWithinBoundingBox(rPoint);
    }

    std::pair<bool, bool> GeometryQuery::IsInside( const Ray_AABB_primitive& rRay ) const {
        if( mMeshIsClosed ){
            return IsInsideClosed(rRay);
        } else {
            return IsInsideOpen(rRay);
        }
    }

    bool GeometryQuery::DoIntersect(PointView rLowerBound, PointView rUpperBound, double Tolerance ) const {
        AABB_primitive aabb(rLowerBound, rUpperBound);
        auto potential_intersections = mTree.Query(aabb);

        const double snap_tolerance = RelativeSnapTolerance(rLowerBound, rUpperBound, Tolerance);

		bool result = false;
		mTriangleMesh.VisitEachTriangle<WithoutNormals>(potential_intersections, [&aabb, &snap_tolerance, &result](const auto& rTriangle){
				if( aabb.intersect(rTriangle.P1, rTriangle.P2, rTriangle.P3, snap_tolerance) ){
					result = true;
					return TriangleMeshView::VisitToken::stop_loop;
				}
				return TriangleMeshView::VisitToken::continue_loop;
		});
        return result;
    }

    std::vector<IndexType> GeometryQuery::GetIntersectedTriangleIds(
            PointView rLowerBound, PointView rUpperBound, double Tolerance ) const {

        // Perform fast search based on aabb tree. Conservative search.
        AABB_primitive aabb(rLowerBound, rUpperBound);
        auto potential_intersections = mTree.Query(aabb);

        std::vector<IndexType> intersected_triangle_ids;
        intersected_triangle_ids.reserve(potential_intersections.size());

        IndexType local_id = 0;
        mTriangleMesh.VisitEachTriangle<WithoutNormals>(potential_intersections, [&](const auto& tri) {
            const auto triangle_id = potential_intersections[local_id++];
            if( aabb.intersect(tri.P1, tri.P2, tri.P3, Tolerance) ){
                intersected_triangle_ids.push_back(triangle_id);
            }
        });

        return intersected_triangle_ids;
    }

    std::pair<bool, bool> GeometryQuery::IsInsideOpen( const Ray_AABB_primitive& rRay ) const {
        double min_distance = MAXD;
        bool is_inside = false;
        auto potential_intersections = mTree.Query(rRay);
        std::optional<std::pair<bool, bool>> final_result;

        mTriangleMesh.VisitEachTriangle<WithoutNormals>(potential_intersections, [&](const auto& tri) {
            double t, u, v;
            bool back_facing, parallel;
            if( rRay.intersect(tri.P1, tri.P2, tri.P3, t, u, v, back_facing, parallel) ) {
                if( !parallel ){
                    double sum_u_v = u+v;
                    if( t < ZEROTOL ){ // origin lies on boundary
                        final_result = std::make_pair(false, true);
                        return TriangleMeshView::VisitToken::stop_loop;
                    }
                    // Ray shoots through boundary.
                    if( u < 0.0+ZEROTOL || v < 0.0+ZEROTOL || sum_u_v > 1.0-ZEROTOL ){
                        final_result = std::make_pair(false, false);
                        return TriangleMeshView::VisitToken::stop_loop;
                    }
                    if( t < min_distance ){
                        is_inside = back_facing;
                        min_distance = t;
                    }
                }
            }
			return TriangleMeshView::VisitToken::continue_loop;
        });

        if (final_result.has_value()) {
            return *final_result;
        }

        return std::make_pair(is_inside, true);
    }

    std::pair<bool, bool> GeometryQuery::IsInsideClosed( const Ray_AABB_primitive& rRay ) const {
        // Get potential ray intersections from AABB tree.
        auto potential_intersections = mTree.Query(rRay);
        IndexType intersection_count = 0;
        std::optional<std::pair<bool, bool>> final_result;

        mTriangleMesh.VisitEachTriangle<WithoutNormals>(potential_intersections, [&](const auto& tri) {
            double t, u, v;
            bool back_facing, parallel;
            // Test for actual intersection, if area is area is not zero.
            if( TriangleUtilities::Area(tri) > 100.0*ZEROTOL
                    && rRay.intersect(tri.P1, tri.P2, tri.P3, t, u, v, back_facing, parallel) ) {
                intersection_count++;
                const double sum_u_v = u+v;
                if( t < ZEROTOL ){ // Origin lies on boundary
                    final_result = std::make_pair(false, true);
                    return TriangleMeshView::VisitToken::stop_loop;
                }
                // Ray shoots through boundary of triangle.
                if( u < 0.0+ZEROTOL || v < 0.0+ZEROTOL || sum_u_v > 1.0-ZEROTOL ){
                    final_result = std::make_pair(false, false);
                    return TriangleMeshView::VisitToken::stop_loop;
                }
                if( parallel ){ // Triangle is parallel to ray.
                    final_result = std::make_pair(false, false);
                    return TriangleMeshView::VisitToken::stop_loop;
                }
            }
            return TriangleMeshView::VisitToken::continue_loop;
        });

        if (final_result.has_value()) {
            return *final_result;
        }
        if( intersection_count % 2 == 1){ // If intersection count is odd, return true.
            return std::make_pair(true, true);
        } else {
            return std::make_pair(false, true);
        }
    }

} // End queso namespace
