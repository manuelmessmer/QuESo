"""Kratos regression tests for trimmed cantilever quadrature rules."""

from __future__ import annotations

import math
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

KM = pytest.importorskip("KratosMultiphysics")

from pyqueso.kratos_interface import Analysis

pytestmark = pytest.mark.kratos


CASES = [
    pytest.param("QuESoSettings1.json", 0.002, None, id="gauss"),
    pytest.param("QuESoSettings2.json", 0.015, None, id="coarse-grid"),
    pytest.param("QuESoSettings3.json", 0.0005, 2592, id="gauss-point-count"),
    pytest.param("QuESoSettings4.json", 0.0005, None, id="ggq-optimal"),
    pytest.param("QuESoSettings5.json", 0.0005, None, id="ggq-reduced-one"),
    pytest.param("QuESoSettings6.json", 0.0005, None, id="ggq-reduced-two"),
    pytest.param("QuESoSettings7.json", 0.0008, None, id="cubic-gauss"),
    pytest.param("QuESoSettings8.json", 0.0008, None, id="cubic-reduced"),
]


def _check_displacement(
    lower: list[float], upper: list[float], geometry, tolerance: float
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
    assert max(errors) < tolerance


@pytest.mark.parametrize(("settings_name", "tolerance", "expected_points"), CASES)
def test_cantilever_quadrature(
    settings_name: str,
    tolerance: float,
    expected_points: int | None,
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("cantilever/trimmed")
    analysis = Analysis(
        queso_settings_path=directory / settings_name,
        analysis_settings_path=directory / "AnalysisSettings.json",
        kratos_parameters_path=directory / "KratosParameters.json",
    )
    analysis.run()
    model_part = analysis.kratos_model.GetModelPart("NurbsMesh")
    geometry = model_part.GetGeometry("NurbsVolume")
    grid = analysis.queso_model.settings("main")["background_grid_settings"]
    _check_displacement(
        grid["lower_bound_xyz"], grid["upper_bound_xyz"], geometry, tolerance
    )

    if expected_points is not None:
        inside_points = sum(
            len(element.integration_points)
            for element in analysis.queso_model.elements("main")
            if not element.is_trimmed
        )
        assert inside_points == expected_points
