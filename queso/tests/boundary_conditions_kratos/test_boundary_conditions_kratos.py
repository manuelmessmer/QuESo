"""Kratos integration tests for ordinary QuESo condition surfaces."""

import unittest
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from pyqueso.kratos_interface import Analysis

DIRECTORY = Path(__file__).parent
PROJECT_ROOT = Path(__file__).resolve().parents[3]


class TestBoundaryConditionsKratos(unittest.TestCase):
    """Verify conditions assembled through the public Analysis API."""

    def _run(self, settings_name: str):
        analysis = Analysis(
            queso_settings_path=DIRECTORY / settings_name,
            analysis_settings_path=DIRECTORY / "AnalysisSettings.json",
            kratos_parameters_path=DIRECTORY / "KratosParameters.json",
        )
        analysis.run()
        return analysis.kratos_model.GetModelPart("Structure")

    def _surface_area(self, model_part) -> float:
        return sum(
            0.5 * condition.GetGeometry().DeterminantOfJacobian()[0]
            for condition in model_part.Conditions
        )

    def test_penalty_support(self) -> None:
        model_part = self._run("QuESoSettings_Penalty.json")
        self.assertGreater(model_part.NumberOfConditions(), 0)
        self.assertAlmostEqual(self._surface_area(model_part), 1183.54304, places=5)
        for condition in model_part.Conditions:
            self.assertAlmostEqual(
                condition.Properties.GetValue(IGA.PENALTY_FACTOR), 1.0e10
            )
            self.assertEqual(list(condition.GetValue(KM.DISPLACEMENT)), [0.0, 0.0, 1.0])

    def test_lagrange_support(self) -> None:
        model_part = self._run("QuESoSettings_Lagrange.json")
        self.assertGreater(self._surface_area(model_part), 0.0)
        for condition in model_part.Conditions:
            self.assertEqual(list(condition.GetValue(KM.DISPLACEMENT)), [0.0, 0.3, 0.0])

    def test_surface_load(self) -> None:
        model_part = self._run("QuESoSettings_SurfaceLoad.json")
        total_force = [0.0, 0.0, 0.0]
        for condition in model_part.Conditions:
            total_force[0] += condition.GetValue(SMA.POINT_LOAD_X)
            total_force[1] += condition.GetValue(SMA.POINT_LOAD_Y)
            total_force[2] += condition.GetValue(SMA.POINT_LOAD_Z)
        self.assertEqual(len(model_part.Conditions), model_part.NumberOfConditions())
        for force in total_force:
            self.assertAlmostEqual(force, 959.49131, places=5)

    def test_pressure_load(self) -> None:
        model_part = self._run("QuESoSettings_Pressure.json")
        total_force = [0.0, 0.0, 0.0]
        for condition in model_part.Conditions:
            total_force[0] += condition.GetValue(SMA.POINT_LOAD_X)
            total_force[1] += condition.GetValue(SMA.POINT_LOAD_Y)
            total_force[2] += condition.GetValue(SMA.POINT_LOAD_Z)
        self.assertAlmostEqual(total_force[0], 0.0, places=5)
        self.assertAlmostEqual(total_force[1], 0.0, places=5)
        self.assertAlmostEqual(total_force[2], -577.978141 * 2.0, places=5)

    def test_coupling_penalty(self) -> None:
        directory = PROJECT_ROOT / "examples/kratos_analysis/cantilever_coupled"
        analysis = Analysis(
            queso_settings_path=directory / "QuESoSettings.json",
            analysis_settings_path=directory / "AnalysisSettings.json",
            kratos_parameters_path=directory / "KratosParameters.json",
        )
        analysis.run()

        root = analysis.kratos_model.GetModelPart("Structure")
        left = root.GetSubModelPart("left")
        right = root.GetSubModelPart("right")
        coupling_conditions = [
            condition for condition in left.Conditions if "Coupling" in condition.Info()
        ]

        self.assertGreater(len(coupling_conditions), 0)
        self.assertFalse(
            any("Coupling" in condition.Info() for condition in right.Conditions)
        )
        root.GetGeometry("LeftVolume")
        root.GetGeometry("RightVolume")
        for condition in coupling_conditions:
            self.assertAlmostEqual(
                condition.Properties.GetValue(IGA.PENALTY_FACTOR), 1.0e10
            )
            self.assertFalse(condition.Properties.GetValue(IGA.COUPLING_SLIP))


if __name__ == "__main__":
    unittest.main()
