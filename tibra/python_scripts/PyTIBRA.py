import TIBRA_PythonApplication as TIBRA_Application
from tibra.python_scripts.b_spline_volume import BSplineVolume
from tibra.python_scripts.helper import *
import os
import shutil

try:
    # TODO: Move to analysis
    import KratosMultiphysics as KM
    kratos_available = True
except:
    kratos_available = False

if kratos_available:
    from kratos_interface.kratos_analysis import Analysis
    from kratos_interface.model_part_utilities import ModelPartUtilities

class PyTIBRA:
    """Main TIBRA python class.

    Provides interface to run TIBRA.
    """
    def __init__(self, json_filename):
        """The constructor"""
        self.parameters = ReadParameters(json_filename)
        self.b_spline_volume = BSplineVolume(self.parameters)
        if self.parameters.EchoLevel() > 0:
            folder_path = "./output/"
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)
            os.mkdir(folder_path)

        self.tibra = TIBRA_Application.TIBRA(self.parameters)

    def Run(self, kratos_model_part = ""):
        if kratos_available and kratos_model_part != "":
            ModelPartUtilities.CreateTIBRAInput(kratos_model_part, self.parameters)
        self.tibra.Run()
        self.elements = self.tibra.GetElements()
        self.conditions = self.tibra.GetConditions()

    def Clear(self):
        self.tibra.Clear()
        self.elements = ""
        self.conditions = ""

    def GetElements(self):
        return self.elements

    def GetConditions(self):
        return self.conditions

    def GetNumberElements(self):
        return self.parameters.NumberOfElements()

    def GetLowerBound(self):
        return self.parameters.LowerBound()

    def GetUpperBound(self):
        return self.parameters.UpperBound()

    def GetBSplineVolume(self):
        return self.b_spline_volume

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

    def GetAnalysis(self):
        return self.analysis

    def ClosestDistances(self, points, directions):
        return self.tibra.ClosestDistances(points, directions)

    #########################################
    #### Kratos related member functions ####
    #########################################

    def UpdateKratosNurbsVolumeModelPart(self, kratos_model_part):
        if kratos_available:
            ModelPartUtilities.RemoveAllElements(kratos_model_part)
            ModelPartUtilities.RemoveAllConditions(kratos_model_part)
            ModelPartUtilities.AddElementsToModelPart(kratos_model_part, self.elements)
            ModelPartUtilities.AddConditionsToModelPart(kratos_model_part, self.conditions, self.GetLowerBound(), self.GetUpperBound())
        else:
            raise Exception("UpdateKratosNurbsVolumeModelPart :: Kratos is not available.")


    def RunKratosAnalysis(self, kratos_settings="KratosParameters.json"):
        if kratos_available:
            self.analysis = Analysis( self.parameters, kratos_settings, self.elements, self.conditions, self.tibra.GetTriangleMesh())
        else:
            raise Exception("RunKratosAnalysis :: Kratos is not available.")