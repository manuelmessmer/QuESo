"""Numerical regression test for the coupled cantilever."""

import math
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

KM = pytest.importorskip("KratosMultiphysics")

from pyqueso.kratos_interface import Analysis
from pyqueso.scripts.helper import point_from_global_to_param_space

pytestmark = pytest.mark.kratos


def _displacement(coordinate: float) -> float:
    inertia = math.pi / 4.0
    length = 10.0
    young_modulus = 100.0
    poisson_ratio = 0.0
    load = -0.1 * math.pi
    shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio))
    shear_factor = (6.0 + 12.0 * poisson_ratio + 6.0 * poisson_ratio**2) / (
        7.0 + 12.0 * poisson_ratio + 4.0 * poisson_ratio**2
    )
    return -(
        load
        * coordinate**2
        * (3.0 * length - coordinate)
        / (6.0 * young_modulus * inertia)
        + load * coordinate / (shear_modulus * math.pi * shear_factor)
    )


def _displacement_from_volume(
    volume,
    coordinate: float,
    bounds_xyz: tuple[list[float], list[float]],
    bounds_uvw: tuple[list[float], list[float]],
) -> float:
    parameter = KM.Vector(
        list(
            point_from_global_to_param_space(
                (0.0, 0.0, coordinate), bounds_xyz, bounds_uvw
            )
        )
    )
    return volume.GlobalCoordinates(parameter)[1]


def test_displacement_profile_and_interface_continuity(
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("coupled_cantilever")
    analysis = Analysis(
        queso_settings_path=directory / "QuESoSettings.json",
        analysis_settings_path=directory / "AnalysisSettings.json",
        kratos_parameters_path=directory / "KratosParameters.json",
    )
    analysis.run()

    root = analysis.kratos_model.GetModelPart("Structure")
    left_volume = root.GetGeometry("LeftVolume")
    right_volume = root.GetGeometry("RightVolume")
    grid = analysis.queso_model.settings("left")["background_grid_settings"]
    bounds_xyz = (grid["lower_bound_xyz"], grid["upper_bound_xyz"])
    bounds_uvw = (grid["lower_bound_uvw"], grid["upper_bound_uvw"])

    reference = _displacement(10.0)
    errors = []
    for coordinate in np.arange(0.0, 5.001, 0.1):
        displacement = _displacement_from_volume(
            left_volume, coordinate, bounds_xyz, bounds_uvw
        )
        errors.append(abs(displacement - _displacement(coordinate)) / reference)
    for coordinate in np.arange(5.0, 10.001, 0.1):
        displacement = _displacement_from_volume(
            right_volume, coordinate, bounds_xyz, bounds_uvw
        )
        errors.append(abs(displacement - _displacement(coordinate)) / reference)

    left_interface = _displacement_from_volume(left_volume, 5.0, bounds_xyz, bounds_uvw)
    right_interface = _displacement_from_volume(
        right_volume, 5.0, bounds_xyz, bounds_uvw
    )
    assert max(errors) < 0.005
    assert left_interface == pytest.approx(right_interface, rel=0.0, abs=5.0e-9)
