// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CGAL_UTILTIES_H
#define CGAL_UTILTIES_H

//// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
//// Project includes
#include "containers/triangle_mesh.h"

namespace tibra {
namespace cgal {

class CGALUtilities
{
public:
    typedef CGAL::Exact_predicates_inexact_constructions_kernel CGALKernalType;
    typedef CGAL::Surface_mesh<CGALKernalType::Point_3> CGALMeshType;
    typedef CGALKernalType::Point_3 CGALPointType;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    static bool CopyMesh(const CGALMeshType& rInputMesh,
                         TriangleMesh& rOutputMesh );

    static bool CopyMesh(const TriangleMesh& rInputMesh,
                         CGALMeshType& rOutputMesh );
};

} // End namespace cgal
} // End namespace tibra

#endif // End CGAL_UTILTIES_H