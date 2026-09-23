"""Regression tests for generated B-spline volumes."""

import json
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest
from pyqueso import Model

CASES = [
    pytest.param(
        f"QuESoSettings{index}.json", f"results_{index}.json", id=f"case-{index}"
    )
    for index in range(1, 9)
]


def _assert_open_and_non_open_splines_match(open_spline, non_open_spline) -> None:
    open_knots = open_spline.t
    non_open_knots = non_open_spline.t
    expected_spacing = non_open_knots[1] - non_open_knots[0]
    for index in range(len(non_open_knots) - 1):
        assert non_open_knots[index + 1] - non_open_knots[index] == pytest.approx(
            expected_spacing, rel=0.0, abs=1.0e-12
        )

    for parameter in np.arange(open_knots[0], open_knots[-1], 100):
        np.testing.assert_allclose(
            non_open_spline(parameter), open_spline(parameter), rtol=0.0, atol=1.0e-12
        )

    for open_knot in open_knots:
        assert any(np.allclose(open_knot, knot) for knot in non_open_knots)


@pytest.mark.parametrize(("settings_name", "results_name"), CASES)
def test_b_spline_volume(
    settings_name: str,
    results_name: str,
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("b_spline_volume")
    model = Model(directory / settings_name)
    volume_open = model.b_spline_volume("main", "open_knot_vector")
    control_points = volume_open.control_points
    knots_u = volume_open.knots_u
    knots_v = volume_open.knots_v
    knots_w = volume_open.knots_w
    polynomial_order = volume_open.polynomial_order
    results = json.loads((directory / results_name).read_text(encoding="utf-8"))

    num_u = volume_open.num_control_points_u
    num_v = volume_open.num_control_points_v
    num_w = volume_open.num_control_points_w
    num_control_points = num_u * num_v * num_w
    assert num_control_points == len(results["cps"])
    assert len(knots_u) == num_u + polynomial_order[0] + 1
    assert len(knots_v) == num_v + polynomial_order[1] + 1
    assert len(knots_w) == num_w + polynomial_order[2] + 1

    assert len(control_points) == num_control_points
    np.testing.assert_allclose(
        np.asarray(control_points), np.asarray(results["cps"]), rtol=0.0, atol=1.0e-12
    )

    control_point_matrix = volume_open.control_points_matrix
    assert control_point_matrix.shape[:3] == (num_u, num_v, num_w)
    count = 0
    for index_w in range(num_w):
        for index_v in range(num_v):
            for index_u in range(num_u):
                np.testing.assert_allclose(
                    control_point_matrix[index_u, index_v, index_w],
                    control_points[count],
                    rtol=0.0,
                    atol=1.0e-12,
                )
                count += 1

    np.testing.assert_allclose(knots_u, results["knots_u"], rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(knots_v, results["knots_v"], rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(knots_w, results["knots_w"], rtol=0.0, atol=1.0e-12)

    volume_non_open = model.b_spline_volume("main", "non_open_knot_vector")
    for direction in range(3):
        _assert_open_and_non_open_splines_match(
            volume_open.spline(direction), volume_non_open.spline(direction)
        )
