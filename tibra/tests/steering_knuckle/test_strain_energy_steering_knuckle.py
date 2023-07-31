# Project imports
from TIBRA_PythonApplication.PyTIBRA import PyTIBRA
import KratosMultiphysics
import unittest

class TestStrainEnergySteeringKnuckleKratos(unittest.TestCase):
    def test_1(self):
        pytibra = PyTIBRA("tibra/tests/steering_knuckle/TIBRAParameters.json")
        pytibra.Run()

        pytibra.RunKratosAnalysis("tibra/tests/steering_knuckle/KratosParameters.json")

        analysis = pytibra.GetAnalysis()
        model_part = analysis.GetModelPart()
        strain_energy = 0.0
        for element in model_part.Elements:
            values = element.CalculateOnIntegrationPoints(KratosMultiphysics.STRAIN_ENERGY, model_part.ProcessInfo)
            weights = element.CalculateOnIntegrationPoints(KratosMultiphysics.INTEGRATION_WEIGHT, model_part.ProcessInfo)
            for value, weight in zip(values, weights):
                strain_energy += value*weight
        self.assertAlmostEqual(strain_energy, 22.67955818930111, 10)

if __name__ == "__main__":
    unittest.main()



