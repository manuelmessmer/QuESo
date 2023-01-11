// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIMMED_DOMAIN_INCLUDE_H
#define TRIMMED_DOMAIN_INCLUDE_H

/// STL includes
#include <memory>
/// Project includes
#include "embedding/trimmed_domain_base.h"
#include "embedding/aabb_tree.h"
#include "utilities/mesh_utilities.h"
#include "embedding/trimmed_domain_on_plane.h"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  TrimmedDomain
 * @author Manuel Messmer
 * @brief  Provides geometrical operations for clipped Brep models (clipped triangle meshes).
 * @details Uses AABB Tree for fast search.
*/
class TrimmedDomain : public TrimmedDomainBase {

public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@brief Builds AABB tree for given mesh.
    ///@param pClippedTriangleMesh
    ///@note mpTriangleMesh must be passed to mTree() and not mpTriangleMesh(), since ptr is moved!
    TrimmedDomain(TriangleMeshPtrType pTriangleMesh, const PointType& rLowerBound, const PointType& rUpperBound, const Parameters& rParameters )
        : TrimmedDomainBase(std::move(pTriangleMesh), rLowerBound, rUpperBound, rParameters), mTree(GetTriangleMesh())
    {
        ///TODO: Improve this!
        const auto& mesh = GetTriangleMesh();
        mClippedMesh.Reserve(mesh.NumOfTriangles());
        MeshUtilities::Append(mClippedMesh, mesh);

        // Construct trimmed domain on plane upper bound of AABB.
        bool upper_bound = true;
        auto p_trimmed_domain_upper_x = MakeUnique<TrimmedDomainOnPlane>(0, upper_bound, mLowerBound, mUpperBound, this);
        auto p_trimmed_domain_upper_y = MakeUnique<TrimmedDomainOnPlane>(1, upper_bound, mLowerBound, mUpperBound, this);
        auto p_trimmed_domain_upper_z = MakeUnique<TrimmedDomainOnPlane>(2, upper_bound, mLowerBound, mUpperBound, this);
        // Construct trimmed domain on plane lower bound of AABB.
        upper_bound = false;
        auto p_trimmed_domain_lower_x = MakeUnique<TrimmedDomainOnPlane>(0, upper_bound, mLowerBound, mUpperBound, this);
        auto p_trimmed_domain_lower_y = MakeUnique<TrimmedDomainOnPlane>(1, upper_bound, mLowerBound, mUpperBound, this);
        auto p_trimmed_domain_lower_z = MakeUnique<TrimmedDomainOnPlane>(2, upper_bound, mLowerBound, mUpperBound, this);

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
        }
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TrimmedDomain. Expects point to be inside AABB. Check is omitted.
    ///@brief Performs ray tracing in direction of the first triangle. Search for all intersection of ray. Inside/Outside is detected
    ///       based on the orientation of the closest intersected triangle (forward or backward facing).
    ///@param rPoint
    ///@return bool
    bool IsInsideTrimmedDomain(const PointType& rPoint) const override;

    ///@brief Triangulates trimmed domain (Surface mesh of outer hull) and return boundary integration points.
    ///@return BoundaryIPVectorPtrType. Boundary integration points to be used for ConstantTerms::Compute.
    BoundaryIPVectorPtrType pGetBoundaryIps() const;

    /// @brief Returns bounding box of trimmed domain. (Might be smaller than the actual domain of element.)
    /// @return BoundingBox (std::pair: first - lower_bound, second - upper_bound)
    const BoundingBox GetBoundingBoxOfTrimmedDomain() const override;

    ///@}
private:

    ///@}
    ///@name Private Members
    ///@{

    AABB_tree mTree;
    TriangleMesh mClippedMesh;
    ///@}
};

///@}

} // End namespace tibra

#endif // TRIMMED_DOMAIN_INCLUDE_H