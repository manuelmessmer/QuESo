import QuESo_PythonApplication as QuESo_App
from queso.python_scripts.b_spline_volume import BSplineVolume
from queso.python_scripts.helper import *
from queso.python_scripts.json_io import JsonIO
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

class PyQuESo:
    """Main QuESo python class.

    Provides interface to run QuESo.
    """
    def __init__(self, json_filename):
        """The constructor"""
        self.settings = JsonIO.ReadSettings(json_filename)
        write_output_to_file = self.settings["general_settings"].GetBool("write_output_to_file")
        output_directory_name = self.settings["general_settings"].GetString("output_directory_name")
        if write_output_to_file:
            folder_path = "./" + output_directory_name + '/'
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)
            os.mkdir(folder_path)

    def Run(self):
        """Run QuESo"""
        self.embedded_model = QuESo_App.EmbeddedModel(self.settings)
        self.embedded_model.CreateAllFromSettings()

    def GetElements(self):
        """ Returns active elements"""
        return self.embedded_model.GetElements()

    def GetConditions(self):
        """ Returns conditions
        """
        return self.embedded_model.GetConditions()

    def GetSettings(self):
        """ Returns settings dictionary
        """
        return self.embedded_model.GetSettings()

    def GetModelInfo(self):
        """ Return model info dictionary
        """
        return self.embedded_model.GetModelInfo()

    def GetBSplineVolume(self, knot_vector_type):
        return BSplineVolume(self.settings, knot_vector_type)

    def GetIntegrationPoints(self):
        integration_points = QuESo_App.IntegrationPointVector()
        # Gather all poitnts (TODO: make this in C++)
        for element in self.GetElements():
            if element.IsTrimmed():
                for point_trimmed_reduced in element.GetIntegrationPoints():
                    weight = point_trimmed_reduced.Weight()
                    if( weight > 0):
                        integration_points.append(point_trimmed_reduced)
            else:
                for point_inside in element.GetIntegrationPoints():
                    integration_points.append(point_inside)
        return integration_points

    def GetAnalysis(self):
        return self.analysis

    #########################################
    #### Kratos related member functions ####
    #########################################
    def UpdateKratosNurbsVolumeModelPart(self, kratos_model_part):
        if kratos_available:
            ModelPartUtilities.RemoveAllElements(kratos_model_part)
            ModelPartUtilities.RemoveAllConditions(kratos_model_part)
            ModelPartUtilities.AddElementsToModelPart(kratos_model_part, self.GetElements())
            ModelPartUtilities.AddConditionsToModelPart(kratos_model_part, self.conditions, self.GetBoundsXYZ(), self.GetBoundsUVW())
        else:
            raise Exception("UpdateKratosNurbsVolumeModelPart :: Kratos is not available.")


    def RunKratosAnalysis(self, kratos_settings="KratosParameters.json"):
        if kratos_available:
            self.analysis = Analysis( self.settings, kratos_settings, self.GetElements(), self.GetConditions() )
        else:
            raise Exception("RunKratosAnalysis :: Kratos is not available.")
