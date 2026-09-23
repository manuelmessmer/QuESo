"""Kratos strain-energy regression tests for the steering knuckle."""

import unittest
from pathlib import Path

import KratosMultiphysics as KM
from pyqueso.kratos_interface import Analysis

DIRECTORY = Path(__file__).parent


class TestStrainEnergySteeringKnuckleKratos(unittest.TestCase):
    """Verify strain energy for each steering-knuckle quadrature setup."""

    def _run_test(
        self, settings_name: str, expected_energy: float, tolerance: float
    ) -> None:
        analysis = Analysis(
            queso_settings_path=DIRECTORY / settings_name,
            analysis_settings_path=DIRECTORY / "AnalysisSettings.json",
            kratos_parameters_path=DIRECTORY / "KratosParameters.json",
        )
        analysis.run()
        model_part = analysis.kratos_model.GetModelPart("NurbsMesh")
        strain_energy = 0.0
        for element in model_part.Elements:
            values = element.CalculateOnIntegrationPoints(
                KM.STRAIN_ENERGY, model_part.ProcessInfo
            )
            weights = element.CalculateOnIntegrationPoints(
                KM.INTEGRATION_WEIGHT, model_part.ProcessInfo
            )
            strain_energy += sum(
                value * weight for value, weight in zip(values, weights)
            )
        self.assertAlmostEqual(strain_energy, expected_energy, delta=tolerance)

    def test_gauss(self) -> None:
        self._run_test("QuESoSettings1.json", 21.7650178, 0.005)

    def test_ggq_optimal(self) -> None:
        self._run_test("QuESoSettings2.json", 21.7650178, 0.005)

    def test_ggq_reduced_one(self) -> None:
        self._run_test("QuESoSettings3.json", 21.7637366, 0.005)

    def test_ggq_reduced_two(self) -> None:
        self._run_test("QuESoSettings4.json", 21.7650821, 0.005)


if __name__ == "__main__":
    unittest.main()
