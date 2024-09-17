# Project imports
from QuESo_PythonApplication.PyQuESo import PyQuESo
from queso.python_scripts.QuESoUnittest import QuESoTestCase
# Kratos imports
import KratosMultiphysics
# External imports
import unittest

class TestStrainEnergySteeringKnuckleKratos(QuESoTestCase):
    def run_test(self, filename, tolerance):
        pyqueso = PyQuESo(filename)
        pyqueso.Run()
        pyqueso.RunKratosAnalysis("queso/tests/steering_knuckle_kratos/KratosParameters.json")

        analysis = pyqueso.GetAnalysis()
        model_part = analysis.GetModelPart()
        strain_energy = 0.0

        for element in model_part.Elements:
            values = element.CalculateOnIntegrationPoints(KratosMultiphysics.STRAIN_ENERGY, model_part.ProcessInfo)
            weights = element.CalculateOnIntegrationPoints(KratosMultiphysics.INTEGRATION_WEIGHT, model_part.ProcessInfo)
            for value, weight in zip(values, weights):
                strain_energy += value*weight
        self.assertAlmostEqual(strain_energy, 21.787, delta=tolerance)

    def test_1(self):
        self.run_test("queso/tests/steering_knuckle_kratos/QuESoSettings1.json", 0.005)

    def test_2(self):
        self.run_test("queso/tests/steering_knuckle_kratos/QuESoSettings2.json", 0.005)

    def test_3(self):
        self.run_test("queso/tests/steering_knuckle_kratos/QuESoSettings3.json", 0.005)

    def test_4(self):
        self.run_test("queso/tests/steering_knuckle_kratos/QuESoSettings4.json", 0.005)

if __name__ == "__main__":
    unittest.main()



