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
    def __init__(self, json_filename):

        with open(json_filename, 'r') as file:
            self.settings = json.load(file)

        general_settings = self.settings["general_settings"]
        echo_level = general_settings["echo_level"]
        if echo_level > 0:
            folder_path = "./output/"
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)
            os.mkdir(folder_path)



    def Run(self):
        general_settings = self.settings["general_settings"]
        input_filename = general_settings["input_filename"]
        self.post_filename = general_settings["postprocess_filename"]
        echo_level = general_settings["echo_level"]
        embedding_flag = general_settings["embedding_flag"]

        self.mesh_settings = self.settings["mesh_settings"]
        self.lower_point = self.mesh_settings["lower_point"]
        self.upper_point = self.mesh_settings["upper_point"]
        polynomial_order = self.mesh_settings["polynomial_order"]
        self.number_of_knot_spans = self.mesh_settings["number_of_knot_spans"]

        trimmed_quadrature_rule_settings = self.settings["trimmed_quadrature_rule_settings"]
        initial_triangle_edge_length = trimmed_quadrature_rule_settings["initial_triangle_edge_length"]
        min_num_boundary_triangles = trimmed_quadrature_rule_settings["min_num_boundary_triangles"]
        moment_fitting_residual = trimmed_quadrature_rule_settings["moment_fitting_residual"]
        init_point_distribution_factor = trimmed_quadrature_rule_settings["init_point_distribution_factor"]

        non_trimmed_quadrature_rule_settings = self.settings["non_trimmed_quadrature_rule_settings"]
        integration_method = non_trimmed_quadrature_rule_settings["integration_method"]

        self.tibra = TIBRA_Application.TIBRA(input_filename,
                                             self.lower_point, self.upper_point, self.number_of_knot_spans, polynomial_order,
                                             initial_triangle_edge_length,
                                             min_num_boundary_triangles,
                                             moment_fitting_residual,
                                             init_point_distribution_factor,
                                             integration_method,
                                             echo_level,
                                             embedding_flag)
        self.elements = self.tibra.GetElements()

    def GetElements(self):
        return self.elements

    def GetNumberElements(self):
        return self.number_of_knot_spans

    def GetLowerPoint(self):
        return self.lower_point

    def GetUpperPoint(self):
        return self.upper_point

    def GetIntegrationPoints(self):
        integration_points = TIBRA_Application.VectorOfIntegrationPoints()
        # Gather all poitnts (TODO: make this in C++)
        for element in self.elements:
            if element.IsTrimmed():
                for point_trimmed_reduced in element.GetIntegrationPointsTrimmed():
                    weight = point_trimmed_reduced.GetWeight()
                    if( weight > 0):
                        integration_points.append(point_trimmed_reduced)
            else:
                for point_inside in element.GetIntegrationPointsInside():
                    integration_points.append(point_inside)
        return integration_points

    def GetTrianglesOnDirichletBoundary(self, dirichlet_condition):
        dirichlet_triangles = TIBRA_Application.VectorOfTriangles()
        #Loop over all elements
        for element in self.elements:
            if element.IsTrimmed():
                for triangle in element.GetDirichletTriangles(dirichlet_condition):
                    dirichlet_triangles.append(triangle)
        return dirichlet_triangles

    def GetTrianglesOnNeumannBoundary(self, neumann_condition):
        neumann_triangles = TIBRA_Application.VectorOfTriangles()
        #Loop over all elements
        for element in self.elements:
            if element.IsTrimmed():
                print(neumann_condition)
                for triangle in element.GetNeumannTriangles(neumann_condition):
                    neumann_triangles.append(triangle)
        return neumann_triangles

    def RunKratosAnalysis(self, dirichlet_settings, neumann_settings, kratos_settings="KratosParameters.json"):
        if kratos_available:
            integration_points = self.GetIntegrationPoints()
            boundary_conditions = []
            for bc in dirichlet_settings:
                dirichlet_triangles = self.GetTrianglesOnDirichletBoundary(bc[0])
                boundary_conditions.append( PenaltySupport(dirichlet_triangles, self.lower_point, self.upper_point, bc[1]) )
            for bc in neumann_settings:
                neumann_triangles = self.GetTrianglesOnDirichletBoundary(bc[0])
                boundary_conditions.append( SurfaceLoad(neumann_triangles, self.lower_point, self.upper_point,bc[1], False) )

            self.analysis = Analysis( self.mesh_settings, kratos_settings, integration_points, boundary_conditions)

    def PostProcess(self):
        if kratos_available:
            model_part = self.analysis.GetModelPart()
            nurbs_volume = model_part.GetGeometry("NurbsVolume")

            # Mesh points are stored in one consecutive array
            self.tibra.ReadWritePostMesh(self.post_filename)
            raw_mesh_points = self.tibra.GetPostMeshPointsRaw()
            num_points = len(raw_mesh_points)//3
            #print(num_points)
            displacements = []
            for i in range(num_points):
                global_point = [raw_mesh_points[3*i], raw_mesh_points[3*i+1], raw_mesh_points[3*i+2]]
                local_point = FromGlobalToParamSpace(global_point, self.lower_point, self.upper_point)

                local_point_kratos = KM.Vector(3)
                local_point_kratos[0] = local_point[0]
                local_point_kratos[1] = local_point[1]
                local_point_kratos[2] = local_point[2]
                # Evaluate deformed nurbs_volume
                deformed_pos_kratos = nurbs_volume.GlobalCoordinates(local_point_kratos)
                deformed_pos = [0, 0, 0]
                deformed_pos[0] = global_point[0] - deformed_pos_kratos[0]
                deformed_pos[1] = global_point[1] - deformed_pos_kratos[1]
                deformed_pos[2] = global_point[2] - deformed_pos_kratos[2]
                displacements.append( deformed_pos )

            TIBRA_Application.WriteDisplacementToVTK(displacements, "output/results.vtk", True)
            #print(displacements[0])

    def GetAnalysis(self):
        return self.analysis
