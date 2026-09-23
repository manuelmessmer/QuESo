"""Kratos strain-energy regression tests for the steering knuckle."""

from collections.abc import Callable
from pathlib import Path

import pytest

KM = pytest.importorskip("KratosMultiphysics")

from pyqueso.kratos_interface import Analysis

pytestmark = pytest.mark.kratos


@pytest.mark.parametrize(
    ("settings_name", "expected_energy"),
    [
        pytest.param("QuESoSettings1.json", 21.7650178, id="gauss"),
        pytest.param("QuESoSettings2.json", 21.7650178, id="ggq-optimal"),
        pytest.param("QuESoSettings3.json", 21.7637366, id="ggq-reduced-one"),
        pytest.param("QuESoSettings4.json", 21.7650821, id="ggq-reduced-two"),
    ],
)
def test_strain_energy(
    settings_name: str,
    expected_energy: float,
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("steering_knuckle") / "kratos"
    analysis = Analysis(
        queso_settings_path=directory / settings_name,
        analysis_settings_path=directory / "AnalysisSettings.json",
        kratos_parameters_path=directory / "KratosParameters.json",
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
        strain_energy += sum(value * weight for value, weight in zip(values, weights))
    assert strain_energy == pytest.approx(expected_energy, rel=0.0, abs=0.005)
