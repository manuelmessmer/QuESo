import QuESo_PythonApplication as QuESo_Application
from queso.python_scripts.b_spline_volume import BSplineVolume
from queso.python_scripts.helper import *
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
        self.parameters = ReadParameters(json_filename)
        output_directory_name = self.parameters.GetGlobalSettings().GetString("output_directory_name")
        if self.parameters.EchoLevel() > 0:
            folder_path = "./" + output_directory_name + '/'
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)
            os.mkdir(folder_path)

    def Run(self, kratos_model_part = ""):
        """Run QuESo"""
        self.queso = QuESo_Application.QuESo(self.parameters)
        self.b_spline_volume = BSplineVolume(self.parameters)
        if kratos_model_part != "":
            if (not kratos_available):
                raise Exception("PyQuESo :: Trying to read from kratos modelpart, but Kratos is not available.")
            self.__ReadInputFromModelPart(kratos_model_part)
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

    def GetNumberElements(self):
        return self.parameters.NumberOfElements()

    def GetLowerBoundDomainXYZ(self):
        return self.parameters.LowerBoundXYZ()

    def GetUpperBoundDomainXYZ(self):
        return self.parameters.UpperBoundXYZ()

    def GetBoundsXYZ(self):
        return [self.parameters.LowerBoundXYZ(), self.parameters.UpperBoundXYZ()]

    def GetBoundsUVW(self):
        return [self.parameters.LowerBoundUVW(), self.parameters.UpperBoundUVW()]

    def GetBSplineVolume(self):
        return BSplineVolume(self.parameters)

    def GetIntegrationPoints(self):
        integration_points = QuESo_Application.IntegrationPointVector()
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
            self.analysis = Analysis( self.parameters, kratos_settings, self.elements, self.conditions, self.queso.GetTriangleMesh())
        else:
            raise Exception("RunKratosAnalysis :: Kratos is not available.")

    def __ReadInputFromModelPart(self, KratosEmbeddedModelPart):
        # Read main mesh
        global_settings = self.parameters.GetGlobalSettings()
        if global_settings.GetString("input_type") == "kratos_modelpart":
            triangle_mesh = self.queso.GetTriangleMesh()
            model_part_name = global_settings.GetString("input_kratos_modelpart_name")
            model_part = KratosEmbeddedModelPart.GetSubModelPart(model_part_name)
            triangle_mesh.Reserve(model_part.NumberOfElements())
            ModelPartUtilities.ReadTriangleMeshFromModelPart(triangle_mesh, model_part, type="Elements")

        # Read condition related meshes
        for condition_settings in self.parameters.GetConditionsSettingsVector():
            if condition_settings.GetString("input_type") == "kratos_modelpart":
                condition = self.queso.CreateNewCondition(condition_settings)
                triangle_mesh = condition.GetTriangleMesh()
                model_part_name = condition_settings.GetString("input_kratos_modelpart_name")
                triangle_mesh.Reserve(model_part.NumberOfConditions())
                model_part = KratosEmbeddedModelPart.GetSubModelPart(model_part_name)
                ModelPartUtilities.ReadTriangleMeshFromModelPart(triangle_mesh, model_part, type="Conditions")