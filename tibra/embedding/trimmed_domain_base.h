// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIMMED_DOMAIN_BASE_INCLUDE_H
#define TRIMMED_DOMAIN_BASE_INCLUDE_H

/// External includes
#include <memory>

/// Project includes
#include "geometries/boundary_integration_point.h"
#include "geometries/triangle_mesh.h"
#include "utilities/parameters.h"

///@name TIBRA Classes
///@{

/// Forward declarations
class Element;

/**
 * @class  TrimmedDomainBase
 * @author Manuel Messmer
 * @brief  Base class for TrimmedDomain. Stores boundary mesh of trimmed domain.
*/
class TrimmedDomainBase {

public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef std::array<double,3> PointType;
    typedef std::vector<BoundaryIntegrationPoint> BoundaryIPVectorType;
    typedef std::unique_ptr<BoundaryIPVectorType> BoundaryIPVectorPtrType;
    typedef std::unique_ptr<TriangleMesh> TriangleMeshPtrType;
    typedef std::pair<PointType, PointType> BoundingBox;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    TrimmedDomainBase(const PointType& rLowerBound, const PointType& rUpperBound)
        : mLowerBound(rLowerBound), mUpperBound(rUpperBound)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TrimmedDomain.
    ///@param rPoint
    ///@return bool
    virtual bool IsInsideTrimmedDomain(const PointType& rPoint) const = 0;

    ///@brief Returns boundary integration points of TrimmedDomain.
    ///@return BoundaryIPVectorPtrType. Boundary integration points to be used for ConstantTerms::Compute.
    virtual BoundaryIPVectorPtrType GetBoundaryIps() const = 0;

    /// @brief Returns bounding box of trimmed domain. (Might be smaller than the actual domain of element.)
    /// @return BoundingBox (std::pair: first - lower_bound, second - upper_bound)
    virtual const BoundingBox GetBoundingBoxOfTrimmedDomain() const = 0;

    /// @brief Return ptr to Triangle mesh (Raw Ptr)
    /// @return const TriangleMesh*
    const TriangleMesh* const pGetTriangleMesh() const{
        return mpTriangleMesh.get();
    }

    /// @brief Return reference to triangle mesh
    /// @return const TriangleMesh&
    const TriangleMesh& GetTriangleMesh() const{
        return *(mpTriangleMesh.get());
    }

    /// @brief Returns part of triangle mesh that IsInDomain.
    /// @param IsInDomain std::function
    /// @return TriangleMeshPtrType (std::unique_ptr)
    TriangleMeshPtrType pGetTriangleMesh(std::function<bool(double, double,double)> &IsInDomain) const {
        // Get Ids of all triangles that are inside given domain.
        std::vector<IndexType> triangle_ids;
        const IndexType num_triangles = mpTriangleMesh->NumOfTriangles();
        for( IndexType triangle_id = 0; triangle_id < num_triangles; ++triangle_id ){
            const auto& p1 = mpTriangleMesh->P1(triangle_id);
            const auto& p2 = mpTriangleMesh->P2(triangle_id);
            const auto& p3 = mpTriangleMesh->P3(triangle_id);
            if( IsInDomain(p1[0], p1[1], p1[2]) ){
                if( IsInDomain(p2[0], p2[1], p2[2]) ){
                    if( IsInDomain(p3[0], p3[1], p3[2]) ){
                        triangle_ids.push_back(triangle_id);
                    }
                }
            }
        }
        // Copy all triangles in (triangle_ids) to new mesh.
        auto p_new_mesh = std::make_unique<TriangleMesh>();
        p_new_mesh->Append(triangle_ids, *mpTriangleMesh);
        return std::move(p_new_mesh);
    }

    ///@}

protected:
    ///@name Protected member variables
    ///@{
    PointType mLowerBound;
    PointType mUpperBound;
    TriangleMeshPtrType mpTriangleMesh;
    ///@}
};
///@}

#endif // TRIMMED_DOMAIN_BASE_INCLUDE_H