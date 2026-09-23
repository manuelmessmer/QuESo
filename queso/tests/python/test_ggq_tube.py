"""Quadrature regression tests for the embedded tube."""

import json
from collections.abc import Callable
from pathlib import Path

import pytest
from pyqueso import Model

CASES = [
    pytest.param("QuESoSettings1.json", "result_ips_gauss.json", id="gauss"),
    pytest.param(
        "QuESoSettings2.json",
        "result_ips_gauss_reduced1.json",
        id="gauss-reduced-one",
    ),
    pytest.param(
        "QuESoSettings3.json",
        "result_ips_gauss_reduced2.json",
        id="gauss-reduced-two",
    ),
    pytest.param(
        "QuESoSettings4.json", "result_ips_reduced_exact.json", id="ggq-optimal"
    ),
    pytest.param(
        "QuESoSettings5.json",
        "result_ips_reduced_order1.json",
        id="ggq-reduced-one",
    ),
    pytest.param(
        "QuESoSettings6.json",
        "result_ips_reduced_order2.json",
        id="ggq-reduced-two",
    ),
]


@pytest.mark.parametrize(("settings_name", "results_name"), CASES)
def test_integration_points(
    settings_name: str,
    results_name: str,
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("ggq_tube")
    model = Model(directory / settings_name)
    model.create()

    integration_points = {
        element.id: [
            [point.x, point.y, point.z, point.weight]
            for point in element.integration_points
        ]
        for element in model.elements("main")
        if not element.is_trimmed
    }
    reference = json.loads((directory / results_name).read_text(encoding="utf-8"))

    for element_id, reference_points in reference.items():
        points = integration_points[int(element_id)]
        assert len(reference_points) == len(points)
        for reference_point, point in zip(reference_points, points):
            for reference_value, value in zip(reference_point, point):
                assert abs(reference_value - value) / abs(reference_value) < 1.0e-11
