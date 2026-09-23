"""Kratos regression tests for trimmed cantilever quadrature rules."""

import math
import unittest
from pathlib import Path

import KratosMultiphysics as KM
import numpy as np
from pyqueso.kratos_interface import Analysis

DIRECTORY = Path(__file__).parent


class TestTrimmedCantileverKratos(unittest.TestCase):
    """Verify displacement accuracy and point counts through Analysis."""

    def _run_test(self, settings_name: str, tolerance: float):
        analysis = Analysis(
            queso_settings_path=DIRECTORY / settings_name,
            analysis_settings_path=DIRECTORY / "AnalysisSettings.json",
            kratos_parameters_path=DIRECTORY / "KratosParameters.json",
        )
        analysis.run()
        model_part = analysis.kratos_model.GetModelPart("NurbsMesh")
        geometry = model_part.GetGeometry("NurbsVolume")
        grid = analysis.queso_model.settings("main")["background_grid_settings"]
        self._check_displacement(
            grid["lower_bound_xyz"], grid["upper_bound_xyz"], geometry, tolerance
        )
        return analysis.queso_model

    def test_gauss_rule(self) -> None:
        self._run_test("QuESoSettings1.json", 0.002)

    def test_coarse_grid(self) -> None:
        self._run_test("QuESoSettings2.json", 0.015)

    def test_gauss_point_count(self) -> None:
        model = self._run_test("QuESoSettings3.json", 0.0005)
        inside_points = sum(
            len(element.integration_points)
            for element in model.elements("main")
            if not element.is_trimmed
        )
        self.assertEqual(inside_points, 2592)

    def test_ggq_optimal(self) -> None:
        self._run_test("QuESoSettings4.json", 0.0005)

    def test_ggq_reduced_one(self) -> None:
        self._run_test("QuESoSettings5.json", 0.0005)

    def test_ggq_reduced_two(self) -> None:
        self._run_test("QuESoSettings6.json", 0.0005)

    def test_cubic_gauss(self) -> None:
        self._run_test("QuESoSettings7.json", 0.0008)

    def test_cubic_reduced(self) -> None:
        self._run_test("QuESoSettings8.json", 0.0008)

    def _check_displacement(
        self, lower: list[float], upper: list[float], geometry, tolerance: float
    ) -> None:
        inertia = math.pi / 4.0
        length = 10.0
        young_modulus = 100.0
        poisson_ratio = 0.0
        load = -0.1 * math.pi
        shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio))
        shear_factor = (6.0 + 12.0 * poisson_ratio + 6.0 * poisson_ratio**2) / (
            7.0 + 12.0 * poisson_ratio + 4.0 * poisson_ratio**2
        )
        reference = -(
            load * length**3 / (3.0 * young_modulus * inertia)
            + load * length / (shear_modulus * math.pi * shear_factor)
        )
        errors = []
        for coordinate in np.arange(0.0, length + 0.001, 0.1):
            parameter = KM.Vector(3)
            parameter[0] = (0.0 - lower[0]) / abs(lower[0] - upper[0])
            parameter[1] = (0.0 - lower[1]) / abs(lower[1] - upper[1])
            parameter[2] = (coordinate - lower[2]) / abs(lower[2] - upper[2])
            displacement = geometry.GlobalCoordinates(parameter)[1]
            expected = -(
                load
                * coordinate**2
                * (3.0 * length - coordinate)
                / (6.0 * young_modulus * inertia)
                + load * coordinate / (shear_modulus * math.pi * shear_factor)
            )
            errors.append(abs(displacement - expected) / abs(reference))
        self.assertLess(max(errors), tolerance)


if __name__ == "__main__":
    unittest.main()
