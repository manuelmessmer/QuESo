# Import Kratos
import KratosMultiphysics as KM

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.IgaApplication as IgaApplication
import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory
import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication

# Project imports
from kratos_interface.geometry_modeler import GeometryModeler
from kratos_interface.custom_analysis_stage import CustomAnalysisStage
#External imports
import matplotlib.pyplot as plt
import math
import numpy as np


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName

class Analysis():
    def __init__(self, filename, integration_points_embedder, analysis_parameters,
            boundary_conditions, material_props, postprocess_flag, gid_output_dest, static = True):
        for modeler in analysis_parameters["modelers"]:
            if modeler["modeler_name"].GetString() == "NurbsGeometryModeler":
                self.lower_point = modeler["Parameters"]["lower_point"].GetVector()
                self.upper_point = modeler["Parameters"]["upper_point"].GetVector()
                self.number_of_knot_spans = modeler["Parameters"]["number_of_knot_spans"].GetVector()
                self.polynomial_order = modeler["Parameters"]["polynomial_order"].GetVector()
        self.boundary_conditions = boundary_conditions
        self.gid_output_dest = gid_output_dest
        self.static = static
        self.material_props = material_props
        self.postprocess_flag = postprocess_flag
        model_ref = KM.Model()
        GeometryModeler.create_geometry(model_ref, self.lower_point, self.upper_point, self.number_of_knot_spans, self.polynomial_order)
        mp_ref = model_ref.GetModelPart("NurbsMesh")
        self.volume_ref = mp_ref.GetGeometry(1)
        # Run analysis
        analysis = self.solve(filename, integration_points_embedder, analysis_parameters)
        self.volume = self.model_part.GetGeometry("NurbsVolume")
        # Run postprocessingvolume
        if postprocess_flag:
            self.cantilever = True
            analysis.postprocess_mdpa(filename, self.volume, self.volume_ref)

    def GetGeometry(self):
        return self.volume

    def GetModelPart(self):
        return self.model_part

    def solve(self, filename, integration_points_embedder, analysis_parameters):
        model = KM.Model()

        analysis = CustomAnalysisStage(model, analysis_parameters, filename, self.volume_ref, self.lower_point, self.upper_point,
            self.boundary_conditions, integration_points_embedder, self.material_props, self.postprocess_flag, self.gid_output_dest)

        analysis.Run()

        self.model_part = model.GetModelPart('NurbsMesh')
        return analysis

    def GetModelPart(self):
        return self.model_part
