import QuESo_PythonApplication as QuESo_App
from queso.python_scripts.b_spline_volume import BSplineVolume
from queso.python_scripts.helper import *
from queso.python_scripts.json_import import JsonImport
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
        self.settings = JsonImport.ReadSettings(json_filename)
        write_output_to_file = self.settings[QuESo_App.MainSettings.general_settings].GetBool(QuESo_App.GeneralSettings.write_output_to_file)
        output_directory_name = self.settings[QuESo_App.MainSettings.general_settings].GetString(QuESo_App.GeneralSettings.output_directory_name)
        if write_output_to_file:
            folder_path = "./" + output_directory_name + '/'
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)
            os.mkdir(folder_path)

    def Run(self):
        """Run QuESo"""
        self.queso = QuESo_App.QuESo(self.settings)
        self.queso.Run()
        self.elements = self.queso.GetElements()
        self.conditions = self.queso.GetConditions()

    def Clear(self):
        self.queso.Clear()
        self.elements = ""
        self.conditions = ""

    def GetElements(self):
        return self.elements

    def GetConditions(self):
        return self.conditions

    def GetSettings(self):
        return self.settings

    # def GetNumberElements(self):
    #     return self.settings[].NumberOfElements()

    # def GetLowerBoundDomainXYZ(self):
    #     return self.settings[].LowerBoundXYZ()

    # def GetUpperBoundDomainXYZ(self):
    #     return self.settings[].UpperBoundXYZ()

    # def GetBoundsXYZ(self):
    #     return [self.settings[].LowerBoundXYZ(), self.settings[].UpperBoundXYZ()]

    # def GetBoundsUVW(self):
    #     return [self.settings[].LowerBoundUVW(), self.settings[].UpperBoundUVW()]

    def GetBSplineVolume(self, knot_vector_type):
        return BSplineVolume(self.settings, knot_vector_type)

    def GetIntegrationPoints(self):
        integration_points = QuESo_App.IntegrationPointVector()
        # Gather all poitnts (TODO: make this in C++)
        for element in self.elements:
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
            ModelPartUtilities.AddElementsToModelPart(kratos_model_part, self.elements)
            ModelPartUtilities.AddConditionsToModelPart(kratos_model_part, self.conditions, self.GetBoundsXYZ(), self.GetBoundsUVW())
        else:
            raise Exception("UpdateKratosNurbsVolumeModelPart :: Kratos is not available.")


    def RunKratosAnalysis(self, kratos_settings="KratosParameters.json"):
        if kratos_available:
            self.analysis = Analysis( self.settings, kratos_settings, self.elements, self.conditions, self.queso.GetTriangleMesh())
        else:
            raise Exception("RunKratosAnalysis :: Kratos is not available.")

    # def __ReadInputFromModelPart(self, KratosEmbeddedModelPart):
    #     # Read main mesh
    #     global_settings = self.settings[].GetGlobalSettings()
    #     if global_settings.GetString("input_type") == "kratos_modelpart":
    #         triangle_mesh = self.queso.GetTriangleMesh()
    #         model_part_name = global_settings.GetString("input_kratos_modelpart_name")
    #         model_part = KratosEmbeddedModelPart.GetSubModelPart(model_part_name)
    #         triangle_mesh.Reserve(model_part.NumberOfElements())
    #         ModelPartUtilities.ReadTriangleMeshFromModelPart(triangle_mesh, model_part, type="Elements")

    #     # Read condition related meshes
    #     for condition_settings in self.settings[].GetConditionsSettingsVector():
    #         if condition_settings.GetString("input_type") == "kratos_modelpart":
    #             condition = self.queso.CreateNewCondition(condition_settings)
    #             triangle_mesh = condition.GetTriangleMesh()
    #             model_part_name = condition_settings.GetString("input_kratos_modelpart_name")
    #             triangle_mesh.Reserve(model_part.NumberOfConditions())
    #             model_part = KratosEmbeddedModelPart.GetSubModelPart(model_part_name)
    #             ModelPartUtilities.ReadTriangleMeshFromModelPart(triangle_mesh, model_part, type="Conditions")