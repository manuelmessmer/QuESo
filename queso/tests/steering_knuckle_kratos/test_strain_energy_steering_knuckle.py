# Project imports
from QuESoPythonModule.kratos_interface.kratos_analysis import KratosAnalysis
from QuESoPythonModule.scripts.queso_unit_test import QuESoTestCase
# Kratos imports
import KratosMultiphysics as KM
# External imports
import unittest

class TestStrainEnergySteeringKnuckleKratos(QuESoTestCase):
    def run_test(self, filename, tolerance):
        model = KM.Model()
        analysis = KratosAnalysis(
            model,
            queso_settings_filename=filename,
            analysis_settings_filename="queso/tests/steering_knuckle_kratos/AnalysisSettings.json",
            kratos_settings_filename="queso/tests/steering_knuckle_kratos/KratosParameters.json",
        )
        analysis.Run()

        model_part = analysis.GetModelPart()
        strain_energy = 0.0

        for element in model_part.Elements:
            values = element.CalculateOnIntegrationPoints(KM.STRAIN_ENERGY, model_part.ProcessInfo)
            weights = element.CalculateOnIntegrationPoints(KM.INTEGRATION_WEIGHT, model_part.ProcessInfo)
            for value, weight in zip(values, weights):
                strain_energy += value*weight
        self.assertAlmostEqual(strain_energy, 21.777, delta=tolerance)

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


