"""Kratos regression tests comparing full and GGQ cantilever rules."""

import json
from collections.abc import Callable
from pathlib import Path

import pytest

KM = pytest.importorskip("KratosMultiphysics")

from pyqueso.kratos_interface import Analysis

pytestmark = pytest.mark.kratos


def _run_analysis(
    directory: Path,
    cross_elements: int,
    axial_elements: int,
    integration_method: str,
) -> tuple[float, int]:
    settings = json.loads(
        (directory / "QuESoSettings.json").read_text(encoding="utf-8")
    )
    component = settings["components"][0]
    component["background_grid_settings"]["number_of_elements"] = [
        cross_elements,
        cross_elements,
        axial_elements,
    ]
    component["non_trimmed_quadrature_rule_settings"]["integration_method"] = (
        integration_method
    )
    settings_path = directory / f"settings_{integration_method}.json"
    settings_path.write_text(json.dumps(settings), encoding="utf-8")

    analysis = Analysis(
        queso_settings_path=settings_path,
        analysis_settings_path=directory / "AnalysisSettings.json",
        kratos_parameters_path=directory / "KratosParameters.json",
    )
    analysis.run()
    geometry = analysis.kratos_model.GetModelPart("NurbsMesh").GetGeometry(
        "NurbsVolume"
    )
    displacement = geometry.GlobalCoordinates(KM.Vector([0.5, 0.5, 1.0]))[1] - 1.0
    quadrature_points = sum(
        len(element.integration_points)
        for element in analysis.queso_model.elements("main")
    )
    return displacement, quadrature_points


@pytest.mark.parametrize(
    "axial_elements",
    [
        pytest.param(3, id="three-elements"),
        pytest.param(4, id="four-elements"),
        pytest.param(5, id="five-elements"),
        pytest.param(6, id="six-elements"),
        pytest.param(7, id="seven-elements"),
    ],
)
def test_reduced_rule_matches_full_rule(
    axial_elements: int,
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("cantilever/ggq")
    reduced_displacement, reduced_points = _run_analysis(
        directory, 8, axial_elements, "GGQ_Optimal"
    )
    full_displacement, full_points = _run_analysis(
        directory, 8, axial_elements, "Gauss"
    )
    assert reduced_points < full_points
    assert reduced_displacement == pytest.approx(
        full_displacement, rel=0.0, abs=5.0e-11
    )
