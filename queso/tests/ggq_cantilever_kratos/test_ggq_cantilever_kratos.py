"""Kratos regression tests comparing full and GGQ-reduced cantilever rules."""

import json
import tempfile
import unittest
from pathlib import Path

import KratosMultiphysics as KM
from pyqueso.kratos_interface import Analysis

DIRECTORY = Path(__file__).parent


class TestGGQCantileverKratos(unittest.TestCase):
    """Compare displacement and quadrature counts for full and GGQ rules."""

    def _run_analysis(
        self, cross_elements: int, axial_elements: int, integration_method: str
    ) -> tuple[float, int]:
        settings = json.loads((DIRECTORY / "QuESoSettings.json").read_text())
        component = settings["components"][0]
        grid = component["background_grid_settings"]
        grid["number_of_elements"] = [cross_elements, cross_elements, axial_elements]
        component["non_trimmed_quadrature_rule_settings"]["integration_method"] = (
            integration_method
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", dir=DIRECTORY, delete=False
        ) as file:
            json.dump(settings, file)
            settings_path = Path(file.name)
        try:
            analysis = Analysis(
                queso_settings_path=settings_path,
                analysis_settings_path=DIRECTORY / "AnalysisSettings.json",
                kratos_parameters_path=DIRECTORY / "KratosParameters.json",
            )
            analysis.run()
        finally:
            settings_path.unlink(missing_ok=True)

        model_part = analysis.kratos_model.GetModelPart("NurbsMesh")
        geometry = model_part.GetGeometry("NurbsVolume")
        parameter = KM.Vector([0.5, 0.5, 1.0])
        displacement = geometry.GlobalCoordinates(parameter)[1] - 1.0
        quadrature_points = sum(
            len(element.integration_points)
            for element in analysis.queso_model.elements("main")
        )
        return displacement, quadrature_points

    def _compare_full_and_reduced(self, axial_elements: int) -> None:
        reduced_displacement, reduced_points = self._run_analysis(
            8, axial_elements, "GGQ_Optimal"
        )
        full_displacement, full_points = self._run_analysis(8, axial_elements, "Gauss")
        self.assertLess(reduced_points, full_points)
        self.assertAlmostEqual(reduced_displacement, full_displacement, places=10)

    def test_three_axial_elements(self) -> None:
        self._compare_full_and_reduced(3)

    def test_four_axial_elements(self) -> None:
        self._compare_full_and_reduced(4)

    def test_five_axial_elements(self) -> None:
        self._compare_full_and_reduced(5)

    def test_six_axial_elements(self) -> None:
        self._compare_full_and_reduced(6)

    def test_seven_axial_elements(self) -> None:
        self._compare_full_and_reduced(7)


if __name__ == "__main__":
    unittest.main()
