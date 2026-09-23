"""Kratos integration tests for QuESo condition surfaces."""

from collections.abc import Callable
from pathlib import Path

import pytest

KM = pytest.importorskip("KratosMultiphysics")
IGA = pytest.importorskip("KratosMultiphysics.IgaApplication")
SMA = pytest.importorskip("KratosMultiphysics.StructuralMechanicsApplication")

from pyqueso.kratos_interface import Analysis

pytestmark = pytest.mark.kratos


def _run(directory: Path, settings_name: str):
    analysis = Analysis(
        queso_settings_path=directory / settings_name,
        analysis_settings_path=directory / "AnalysisSettings.json",
        kratos_parameters_path=directory / "KratosParameters.json",
    )
    analysis.run()
    return analysis.kratos_model.GetModelPart("Structure")


def _surface_area(model_part) -> float:
    return sum(
        0.5 * condition.GetGeometry().DeterminantOfJacobian()[0]
        for condition in model_part.Conditions
    )


def _boundary_case(copy_test_data: Callable[[str], Path]) -> Path:
    directory = copy_test_data("boundary_conditions")
    copy_test_data("steering_knuckle")
    return directory


def test_penalty_support(copy_test_data: Callable[[str], Path]) -> None:
    model_part = _run(_boundary_case(copy_test_data), "QuESoSettings_Penalty.json")
    assert model_part.NumberOfConditions() > 0
    assert _surface_area(model_part) == pytest.approx(1183.54304, rel=0.0, abs=5.0e-6)
    for condition in model_part.Conditions:
        assert condition.Properties.GetValue(IGA.PENALTY_FACTOR) == pytest.approx(
            1.0e10, rel=0.0, abs=5.0e-8
        )
        assert list(condition.GetValue(KM.DISPLACEMENT)) == [0.0, 0.0, 1.0]


def test_lagrange_support(copy_test_data: Callable[[str], Path]) -> None:
    model_part = _run(_boundary_case(copy_test_data), "QuESoSettings_Lagrange.json")
    assert _surface_area(model_part) > 0.0
    for condition in model_part.Conditions:
        assert list(condition.GetValue(KM.DISPLACEMENT)) == [0.0, 0.3, 0.0]


def test_surface_load(copy_test_data: Callable[[str], Path]) -> None:
    model_part = _run(_boundary_case(copy_test_data), "QuESoSettings_SurfaceLoad.json")
    total_force = [0.0, 0.0, 0.0]
    for condition in model_part.Conditions:
        total_force[0] += condition.GetValue(SMA.POINT_LOAD_X)
        total_force[1] += condition.GetValue(SMA.POINT_LOAD_Y)
        total_force[2] += condition.GetValue(SMA.POINT_LOAD_Z)
    assert len(model_part.Conditions) == model_part.NumberOfConditions()
    assert total_force == pytest.approx([959.49131] * 3, rel=0.0, abs=5.0e-6)


def test_pressure_load(copy_test_data: Callable[[str], Path]) -> None:
    model_part = _run(_boundary_case(copy_test_data), "QuESoSettings_Pressure.json")
    total_force = [0.0, 0.0, 0.0]
    for condition in model_part.Conditions:
        total_force[0] += condition.GetValue(SMA.POINT_LOAD_X)
        total_force[1] += condition.GetValue(SMA.POINT_LOAD_Y)
        total_force[2] += condition.GetValue(SMA.POINT_LOAD_Z)
    assert total_force == pytest.approx(
        [0.0, 0.0, -577.978141 * 2.0], rel=0.0, abs=5.0e-6
    )


def test_coupling_penalty(copy_test_data: Callable[[str], Path]) -> None:
    directory = copy_test_data("coupled_cantilever")
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

    assert coupling_conditions
    assert not any("Coupling" in condition.Info() for condition in right.Conditions)
    root.GetGeometry("LeftVolume")
    root.GetGeometry("RightVolume")
    for condition in coupling_conditions:
        assert condition.Properties.GetValue(IGA.PENALTY_FACTOR) == pytest.approx(
            1.0e10, rel=0.0, abs=5.0e-8
        )
        if hasattr(IGA, "COUPLING_SLIP"):
            assert not condition.Properties.GetValue(IGA.COUPLING_SLIP)
