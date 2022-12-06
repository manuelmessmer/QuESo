import TIBRA_PythonApplication as TIBRA_Application
from tibra.python_scripts.helper import *
import json
import os
import shutil

try:
    # TODO: Move to analysis
    import KratosMultiphysics as KM
    kratos_available = True
except:
    print("KratosMultiphysics is not available")
    kratos_available = False

if kratos_available:
    from kratos_interface.kratos_analysis import Analysis
    from kratos_interface.weak_bcs import PenaltySupport
    from kratos_interface.weak_bcs import SurfaceLoad
    from kratos_interface.bounding_box_bcs import DirichletCondition
    from kratos_interface.bounding_box_bcs import NeumannCondition

class PyTIBRA:
    """Main TIBRA python class.

    Provides interface to run TIBRA.
    """
    def __init__(self, json_filename):
        """The constructor"""
        self.parameters = ReadParameters(json_filename)

        if self.parameters.EchoLevel() > 0:
            folder_path = "./output/"
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)
            os.mkdir(folder_path)


    def Run(self):
        self.tibra = TIBRA_Application.TIBRA(self.parameters)
        self.elements = self.tibra.GetElements()

    def GetElements(self):
        return self.elements

    def GetNumberElements(self):
        return self.parameters.NumberOfElememnts()

    def GetLowerBound(self):
        return self.parameters.LowerBound()

    def GetUpperBound(self):
        return self.parameters.UpperBound()

    def GetIntegrationPoints(self):
        integration_points = TIBRA_Application.IntegrationPointVector()
        # Gather all poitnts (TODO: make this in C++)
        for element in self.elements:
            if element.IsTrimmed():
                for point_trimmed_reduced in element.GetIntegrationPoints():
                    weight = point_trimmed_reduced.GetWeight()
                    if( weight > 0):
                        integration_points.append(point_trimmed_reduced)
            else:
                for point_inside in element.GetIntegrationPoints():
                    integration_points.append(point_inside)
        return integration_points


    def GetTrianglesOnBoundary(self, condition):
        triangle_mesh = TIBRA_Application.TriangleMesh()
        #Loop over all elements
        for element in self.elements:
            if element.IsTrimmed():
                new_mesh = element.GetBCTriangleMesh(condition)
                triangle_mesh.Append(new_mesh)

        return triangle_mesh

    def RunKratosAnalysis(self, dirichlet_settings, neumann_settings, kratos_settings="KratosParameters.json"):
        if kratos_available:
            integration_points = self.GetIntegrationPoints()
            boundary_conditions = []
            for bc in dirichlet_settings:
                dirichlet_triangles = self.GetTrianglesOnBoundary(bc[0])
                boundary_conditions.append( PenaltySupport(dirichlet_triangles, self.GetLowerBound(), self.GetUpperBound(), bc[1]) )
            for bc in neumann_settings:
                neumann_triangles = self.GetTrianglesOnBoundary(bc[0])
                boundary_conditions.append( SurfaceLoad(neumann_triangles, self.GetLowerBound(), self.GetUpperBound(), bc[1], False) )

            self.analysis = Analysis( self.parameters, kratos_settings, integration_points, boundary_conditions)

    def PostProcess(self):
        if kratos_available:
            model_part = self.analysis.GetModelPart()
            nurbs_volume = model_part.GetGeometry("NurbsVolume")

            # Mesh points are stored in one consecutive array
            self.tibra.ReadWritePostMesh()
            mesh_points = self.tibra.GetPostMeshPoints()

            #print(num_points)
            displacements = TIBRA_Application.PointVector()
            for point in mesh_points:
                global_point = [point[0], point[1], point[2]]
                local_point = FromGlobalToParamSpace(global_point, self.GetLowerBound(), self.GetUpperBound())

                local_point_kratos = KM.Vector(3)
                local_point_kratos[0] = local_point[0]
                local_point_kratos[1] = local_point[1]
                local_point_kratos[2] = local_point[2]
                # Evaluate deformed nurbs_volume
                deformed_pos_kratos = nurbs_volume.GlobalCoordinates(local_point_kratos)
                deformed_pos = TIBRA_Application.Point(deformed_pos_kratos[0] - global_point[0],
                    deformed_pos_kratos[1] - global_point[1],
                    deformed_pos_kratos[2] - global_point[2])

                displacements.append( deformed_pos )

            TIBRA_Application.WriteDisplacementToVTK(displacements, "output/results.vtk", True)

    def GetAnalysis(self):
        return self.analysis
