# Import Kratos
import KratosMultiphysics as KM

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.IgaApplication as IgaApplication
import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory
import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication

# Project imports
from kratos_interface.custom_analysis_stage import CustomAnalysisStage
#External imports
import os

def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName

class Analysis():
    """Wrapper for customized Kratos analysis stage.

    Runs simulation.
    """
    def __init__(self, mesh_settings, kratos_settings_filename, integration_points_embedder, boundary_conditions ):
        """The constructor."""
        self.model = KM.Model()
        analysis = CustomAnalysisStage(self.model, mesh_settings, kratos_settings_filename, integration_points_embedder, boundary_conditions)
        analysis.Run()

    def GetModelPart(self):
        return self.model.GetModelPart("NurbsMesh")
