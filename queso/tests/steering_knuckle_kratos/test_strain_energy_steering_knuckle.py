# Project imports
from TIBRA_PythonApplication.PyTIBRA import PyTIBRA
import KratosMultiphysics
import unittest

class TestStrainEnergySteeringKnuckleKratos(unittest.TestCase):
    def run_test(self, filename, tolerance):
        pytibra = PyTIBRA(filename)
        pytibra.Run()
        pytibra.RunKratosAnalysis("queso/tests/steering_knuckle_kratos/KratosParameters.json")

        analysis = pytibra.GetAnalysis()
        model_part = analysis.GetModelPart()
        strain_energy = 0.0

        for element in model_part.Elements:
            values = element.CalculateOnIntegrationPoints(KratosMultiphysics.STRAIN_ENERGY, model_part.ProcessInfo)
            weights = element.CalculateOnIntegrationPoints(KratosMultiphysics.INTEGRATION_WEIGHT, model_part.ProcessInfo)
            for value, weight in zip(values, weights):
                strain_energy += value*weight
        self.assertAlmostEqual(strain_energy, 21.788415916271443, tolerance)

    def test_1(self):
        self.run_test("queso/tests/steering_knuckle_kratos/TIBRAParameters1.json", 10)

    def test_2(self):
        self.run_test("queso/tests/steering_knuckle_kratos/TIBRAParameters2.json", 10)

    def test_3(self):
        self.run_test("queso/tests/steering_knuckle_kratos/TIBRAParameters3.json", 2)

    def test_4(self):
        self.run_test("queso/tests/steering_knuckle_kratos/TIBRAParameters4.json", 3)

if __name__ == "__main__":
    unittest.main()



