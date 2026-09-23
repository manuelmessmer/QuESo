"""Integration tests for embedded boundary-condition surfaces."""

from collections.abc import Callable
from pathlib import Path

import pyqueso
import pytest

EXPECTED_ACTIVE_AREAS = {
    "SurfaceLoadCondition": 332.37753991575306,
    "PressureLoadCondition": 577.978140756255,
    "LagrangeSupportCondition": 921.1636351534689,
    "PenaltySupportCondition": 1183.5430441558533,
}


def _segmented_area(condition) -> float:
    return sum(
        pyqueso.mesh.area(segment.triangle_mesh) for segment in condition.segments
    )


def _check_values(model, mesh_directory: Path) -> None:
    for condition in model.conditions("main"):
        settings = condition.settings
        condition_type = settings.get_string("condition_type")
        expected_area = EXPECTED_ACTIVE_AREAS[condition_type]

        if condition_type == "SurfaceLoadCondition":
            assert (
                Path(settings.get_string("input_filename")) == mesh_directory / "N1.stl"
            )
            assert settings.get_double("modulus") == pytest.approx(
                5.0, rel=0.0, abs=1.0e-10
            )
            assert settings.get_double_vector("direction") == pytest.approx(
                [-1.0, 2.0, 3.0], rel=0.0, abs=1.0e-10
            )
        elif condition_type == "PressureLoadCondition":
            assert (
                Path(settings.get_string("input_filename")) == mesh_directory / "N2.stl"
            )
            assert settings.get_double("modulus") == pytest.approx(
                2.0, rel=0.0, abs=1.0e-10
            )
        elif condition_type == "LagrangeSupportCondition":
            assert (
                Path(settings.get_string("input_filename")) == mesh_directory / "N3.stl"
            )
            assert settings.get_double_vector("value") == pytest.approx(
                [0.0, 0.3, 0.0], rel=0.0, abs=1.0e-10
            )
        elif condition_type == "PenaltySupportCondition":
            assert (
                Path(settings.get_string("input_filename")) == mesh_directory / "D1.stl"
            )
            assert settings.get_double_vector("value") == pytest.approx(
                [0.0, 0.0, 0.0], rel=0.0, abs=1.0e-10
            )
            assert settings.get_double("penalty_factor") == pytest.approx(
                1.0e10, rel=0.0, abs=1.0e-10
            )
        else:
            pytest.fail(f"Unexpected condition type: {condition_type}")

        assert _segmented_area(condition) == pytest.approx(
            expected_area, rel=0.0, abs=1.0e-5
        )
        assert _segmented_area(condition) == pytest.approx(
            expected_area, rel=0.0, abs=1.0e-5
        )


@pytest.mark.parametrize(
    "settings_name",
    [
        pytest.param("QuESoSettings1.json", id="serial"),
        pytest.param("QuESoSettings2.json", id="partitioned"),
    ],
)
def test_boundary_conditions(
    settings_name: str,
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("boundary_conditions")
    mesh_directory = copy_test_data("steering_knuckle")
    model = pyqueso.Model(directory / settings_name)
    model.create()
    _check_values(model, mesh_directory)


def test_coupling_penalty(copy_test_data: Callable[[str], Path]) -> None:
    directory = copy_test_data("coupled_cantilever")
    model = pyqueso.Model(directory / "QuESoSettings.json")
    model.create()

    left_conditions = {
        condition.settings.get_int("condition_id"): condition
        for condition in model.conditions("left")
    }
    right_condition_ids = {
        condition.settings.get_int("condition_id")
        for condition in model.conditions("right")
    }
    coupling = left_conditions[10]
    settings = coupling.settings

    assert 10 not in right_condition_ids
    assert settings.get_string("condition_type") == "CouplingPenaltyCondition"
    assert settings.get_string("coupling_partner") == "right"
    assert settings.get_double("penalty_factor") == pytest.approx(
        1.0e10, rel=0.0, abs=5.0e-8
    )
    assert not settings.get_bool("slip")
    assert coupling.num_segments > 0
    assert all(segment.is_in_active_element for segment in coupling.segments)

    interface_mesh = pyqueso.mesh.TriangleMesh()
    pyqueso.io.read_mesh_from_stl(interface_mesh, str(directory / "interface.stl"))
    assert _segmented_area(coupling) == pytest.approx(
        pyqueso.mesh.area(interface_mesh), rel=0.0, abs=5.0e-9
    )
