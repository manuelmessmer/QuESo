#ifndef VTK_UTILTIES_H
#define VTK_UTILTIES_H

namespace CGAL{

template<typename PM>
void polygon_mesh_to_vtkUnstructured_(const PM& pmesh,//PolygonMesh
                                      const char* filename);

}

#endif // VTK_UTILTIES_H