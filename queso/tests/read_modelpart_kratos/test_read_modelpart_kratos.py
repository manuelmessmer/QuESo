# Project imports
from QuESo_PythonApplication.PyQuESo import PyQuESo
import KratosMultiphysics as KM
from queso.python_scripts.QuESoUnittest import QuESoTestCase
# External imports
import unittest

class TestReadModelPartKratos(QuESoTestCase):
    def test_kratos_input(self):
        pyqueso = PyQuESo("queso/tests/read_modelpart_kratos/QuESoParameters.json")

        # Read embedded model part
        model = KM.Model()
        embedded_model_part = model.CreateModelPart("hook")
        model_part_io = KM.ModelPartIO("queso/tests/read_modelpart_kratos/input_kratos")
        model_part_io.ReadModelPart(embedded_model_part)

        pyqueso.Run(embedded_model_part)
        pyqueso.RunKratosAnalysis("queso/tests/read_modelpart_kratos/KratosParameters.json")

        analysis = pyqueso.GetAnalysis()
        computing_model_part = analysis.GetModelPart()
        strain_energy = 0.0
        for element in computing_model_part.Elements:
            values = element.CalculateOnIntegrationPoints(KM.STRAIN_ENERGY, computing_model_part.ProcessInfo)
            weights = element.CalculateOnIntegrationPoints(KM.INTEGRATION_WEIGHT, computing_model_part.ProcessInfo)
            for value, weight in zip(values, weights):
                strain_energy += value*weight
        self.assertAlmostEqual(strain_energy, 53588.324, delta=0.01)

if __name__ == "__main__":
    unittest.main()