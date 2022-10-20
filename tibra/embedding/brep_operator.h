// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BREP_OPERATOR_INCLUDE_H
#define BREP_OPERATOR_INCLUDE_H

/// External includes
#include <memory>

/// Project includes
#include "geometries/triangle_mesh.h"
#include "geometries/element.h"
#include "embedding/aabb_tree.h"
#include "embedding/clipper.h"

///@name TIBRA Classes
///@{

/**
 * @class  BRepOperator
 * @author Manuel Messmer
 * @brief  Provides geometrical operations for Brep models.
 * @details Uses AABB Tree for fast search.
*/
class BRepOperator {

public:
    ///@name Type Definitions
    ///@{
    typedef TriangleMesh::Vector3d PointType;

    enum IntersectionStatus {Inside, Outside, Trimmed};
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@brief Builds AABB tree for given mesh.
    ///@param rTriangleMesh
    BRepOperator(const TriangleMesh& rTriangleMesh) : mTriangleMesh(rTriangleMesh), mTree(rTriangleMesh)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TriangleMesh.
    ///@param rPoint
    ///@return bool
    bool IsInside(const PointType& rPoint) const;

    ///@brief Returns intersections state of element.
    ///@param rElement
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    IntersectionStatus GetIntersectionState(const Element& rElement) const;

    ///@brief Returns intersections state of element.
    ///@param rLowerBound
    ///@param rUpperBound
    ///@param Tolerance Tolerance reduces element slightly. If Tolerance=0 touch is detected as intersection.
    ///                 If Tolerance>0, touch is not detected as intersection.
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).

    IntersectionStatus GetIntersectionState(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance) const;



    ///@todo Just return ids. Much more efficient
    std::unique_ptr<TriangleMesh> GetIntersectedTriangles( const TriangleMesh& rTriangleMesh, const PointType& rLowerBound, const PointType& rUpperBound ) const{

        // Perform fast search based on aabb tree. Conservative search.
        AABB_primitive aabb(rLowerBound, rUpperBound);
        auto potential_intersections = mTree.Query(aabb);

        std::vector<IndexType> intersected_triangle_ids{};
        int count = 0;

        for( auto triangle_id : potential_intersections){
            const auto& p1 = rTriangleMesh.P1(triangle_id);
            const auto& p2 = rTriangleMesh.P2(triangle_id);
            const auto& p3 = rTriangleMesh.P3(triangle_id);
            // If tolerance>=0 intersection is not detected.
            const double tolerance_1 = 1e-8;
            // Perform actual intersection test.
            if( aabb.intersect(p1, p2, p3, tolerance_1) ){
                intersected_triangle_ids.push_back(triangle_id);
            }
        }

        // Copy intersected triangles to new mesh.
        TriangleMesh new_mesh{};
        new_mesh.Copy(intersected_triangle_ids, rTriangleMesh);

        return std::make_unique<TriangleMesh>(new_mesh);
    }


    std::unique_ptr<TriangleMesh> ClipTriangleMesh(const PointType& rLowerBound, const PointType& rUpperBound){

        auto p_intersected_triangles = GetIntersectedTriangles(mTriangleMesh, rLowerBound, rUpperBound);

        TriangleMesh new_mesh{};
        std::map<IndexType, IndexType> index_map{};
        IndexType vertex_count = 0;
        for( IndexType triangle_id = 0; triangle_id < p_intersected_triangles->NumOfTriangles(); ++triangle_id ){
            const auto& P1 = p_intersected_triangles->P1(triangle_id);
            const auto& P2 = p_intersected_triangles->P2(triangle_id);
            const auto& P3 = p_intersected_triangles->P3(triangle_id);

            if(    IsContained(P1, rLowerBound, rUpperBound )
                && IsContained(P2, rLowerBound, rUpperBound )
                && IsContained(P3, rLowerBound, rUpperBound ) ){ // Triangle is fully contained, does not need to be clipped.

                new_mesh.AddVertex(P1);
                new_mesh.AddVertex(P2);
                new_mesh.AddVertex(P3);

                // Copy triangles and normals.
                new_mesh.AddTriangle({vertex_count, vertex_count+1, vertex_count+2 });
                vertex_count += 3;
                new_mesh.AddNormal( p_intersected_triangles->Normal(triangle_id) );
            }
            else { // Triangle needs to be clipped.
                //throw std::runtime_error("Wht the fuck!!");
                auto polygon = Clipper::ClipTriangle(P1, P2, P3, rLowerBound, rUpperBound);
                for( auto triangle : polygon->GetTriangles() ){
                    new_mesh.AddVertex(triangle[0]);
                    new_mesh.AddVertex(triangle[1]);
                    new_mesh.AddVertex(triangle[2]);

                    // Copy triangles and normals.
                    new_mesh.AddTriangle({vertex_count, vertex_count+1, vertex_count+2 });
                    vertex_count += 3;
                    new_mesh.AddNormal( p_intersected_triangles->Normal(triangle_id) );
                }
            }

        }
        new_mesh.Check();
        return std::make_unique<TriangleMesh>(new_mesh);
    }

    ///@}


private:
    ///@name Private Operations
    ///@{

    inline bool IsContained(const PointType& rPoint, const PointType& rLowerBound, const PointType& rUpperBound){
        if(    rPoint[0] < rLowerBound[0]
            || rPoint[0] > rUpperBound[0]
            || rPoint[1] < rLowerBound[1]
            || rPoint[1] > rUpperBound[1]
            || rPoint[2] < rLowerBound[2]
            || rPoint[2] > rUpperBound[2] )
        {
            return false;
        }

        return true;
    }

    ///@}
    ///@name Private Members
    ///@{

    AABB_tree mTree;
    const TriangleMesh& mTriangleMesh;
    ///@}
};

///@}

#endif // BREP_OPERATOR_INCLUDE_H