import os
# Import QuESo
import QuESo_PythonApplication as QuESo
# Project imports
from kratos_interface.custom_analysis_stage import CustomAnalysisStage
# Import Kratos
import KratosMultiphysics as KM

def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName

class Analysis():
    """
    Wrapper for customized Kratos analysis stage.

    This class provides a wrapper around the `CustomAnalysisStage` to run the simulation in Kratos.
    It sets up the model, creates the elements, applies the boundary conditions, and executes the analysis.
    """
    def __init__(self,
            settings: QuESo.Settings,
            kratos_settings_filename: str,
            elements: QuESo.ElementVector,
            boundary_conditions: QuESo.ConditionVector
        ) -> None:
        """
        Constructor for the `Analysis` class.

        This constructor initializes the Kratos model and sets up the `CustomAnalysisStage` with the given
        settings, elements, and boundary conditions. It then runs the analysis.

        Args:
            settings (QuESo.Settings): The settings for the analysis, including elements and boundary conditions.
            kratos_settings_filename (str): The path to the Kratos settings file in JSON format.
            elements (QuESo.ElementVector): A list of elements to be added to the model part.
            boundary_conditions (QuESo.ConditionVector): A list of boundary conditions to be applied to the model.
        """
        self.model = KM.Model()
        analysis = CustomAnalysisStage(self.model, settings, kratos_settings_filename, elements, boundary_conditions)
        analysis.Run()

    def GetModelPart(self) -> None:
        """
        Returns the Kratos model part created in the analysis.

        This method retrieves the model part named "NurbsMesh" from the Kratos model.

        Returns:
            ModelPart: The model part named "NurbsMesh" in the Kratos model.
        """
        return self.model.GetModelPart("NurbsMesh")
