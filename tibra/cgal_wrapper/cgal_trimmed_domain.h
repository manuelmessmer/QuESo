// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CGAL_TRIMMED_DOMAIN_INCLUDE_H
#define CGAL_TRIMMED_DOMAIN_INCLUDE_H

/// CGAL includes
// Domain
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>

// Project includes
#include "cgal_wrapper/cgal_utilities.h"
#include "embedding/trimmed_domain_base.h"
#include "containers/triangle_mesh.h"

namespace cgal {

///@name TIBRA Classes
///@{

typedef CGAL::Exact_predicates_inexact_constructions_kernel CGALKernalType;
typedef CGALKernalType::Point_3 CGALPointType;
typedef CGAL::Surface_mesh<CGALPointType> CGALMeshType;

/**
 * @class  CGALTrimmedDomain
 * @author Manuel Messmer
 * @brief TrimmedDomain with CGAL dependency. Stores boundary mesh of trimmed domain.
 *        In constrast to trimmed_domain.h, the boundary mesh is closed.
*/
class CGALTrimmedDomain : public TrimmedDomainBase {

public:
    ///@name Type Definitions
    ///@{

    typedef std::unique_ptr<CGALMeshType> CGALMeshPtrType;
    typedef std::unique_ptr<TriangleMesh> TriangleMeshPtrType;
    typedef CGAL::Side_of_triangle_mesh<CGALMeshType, CGALKernalType> CGALInsideTestType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@param pCGALMesh ptr to closed surface mesh of trimmed domain.
    ///@param rLowerBound of element (full domain, not only the trimmed part).
    ///@param rUpperBound of element (full domain, not only the trimmed part).
    CGALTrimmedDomain(CGALMeshPtrType& pCGALMesh, const PointType& rLowerBound, const PointType& rUpperBound )
        : TrimmedDomainBase(rLowerBound, rUpperBound), mpCGALMesh(std::move(pCGALMesh))
    {
        mpTriangleMesh = std::make_unique<TriangleMesh>();
        cgal::CGALUtilities::CopyMesh(*mpCGALMesh, *mpTriangleMesh);
        mpCGALInsideTest = std::make_unique<CGALInsideTestType>(*mpCGALMesh);
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TrimmedDomain.
    ///@param rPoint
    ///@return bool
    bool IsInsideTrimmedDomain(const PointType& rPoint) const override;

    ///@brief Returns boundary integration points of TrimmedDomain.
    ///@return BoundaryIPVectorPtrType. Boundary integration points to be used for ConstantTerms::Compute.
    BoundaryIPVectorPtrType pGetBoundaryIps() const override;

    /// @brief Returns bounding box of trimmed domain. (Might be smaller than the actual domain of element.)
    /// @return BoundingBox (std::pair: first - lower_bound, second - upper_bound)
    const BoundingBox GetBoundingBoxOfTrimmedDomain() const override;

    ///@}
private:

    ///@name Private member variables
    ///@{
    std::unique_ptr<CGALInsideTestType> mpCGALInsideTest;
    CGALMeshPtrType mpCGALMesh;
    ///@}

}; // End class

} // End namespace cgal

#endif // End