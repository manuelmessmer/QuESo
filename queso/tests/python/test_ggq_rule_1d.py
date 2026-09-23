"""Tests for one-dimensional generalized Gaussian quadrature rules."""

import numpy as np
import pytest
from pyqueso import IntegrationMethod, testing
from scipy.interpolate import BSpline


def _check_ggq_rule(
    points,
    polynomial_degree: int,
    continuity: int,
    num_elements: int,
    must_pass: bool,
) -> None:
    num_dofs = (
        (polynomial_degree + 1) * 2
        + (num_elements - 1) * (polynomial_degree - continuity)
        - polynomial_degree
        - 1
    )
    expected_num_points = int(np.ceil(num_dofs / 2))
    nodes = [point[0] for point in points]
    weights = [point[1] for point in points]

    knots = np.zeros(polynomial_degree + 1)
    knots = np.concatenate(
        (
            knots,
            np.repeat(
                np.linspace(0.0, 1.0, num_elements + 1)[1:-1],
                polynomial_degree - continuity,
            ),
            np.ones(polynomial_degree + 1),
        )
    )
    dimension = len(knots) - polynomial_degree - 1
    target = np.asarray(
        [
            (knots[polynomial_degree + index + 1] - knots[index])
            / (polynomial_degree + 1)
            for index in range(dimension)
        ]
    )
    collocation = BSpline.design_matrix(nodes, knots, polynomial_degree).toarray()
    error = np.linalg.norm((target - np.asarray(weights).dot(collocation)) / dimension)

    if must_pass:
        assert error < 1.0e-15
        assert len(nodes) == expected_num_points
    else:
        assert error > 1.0e-10


CASES = [
    pytest.param(
        2, IntegrationMethod.GGQ_OPTIMAL, 1, (4, 0), [(5, 0, 2)], id="optimal-p2"
    ),
    pytest.param(
        3,
        IntegrationMethod.GGQ_OPTIMAL,
        1,
        (6, 1),
        [(7, 1, 2), (7, 0, 2)],
        id="optimal-p3",
    ),
    pytest.param(
        4,
        IntegrationMethod.GGQ_OPTIMAL,
        1,
        (8, 2),
        [(9, 2, 2), (8, 1, 3)],
        id="optimal-p4",
    ),
    pytest.param(
        2, IntegrationMethod.GGQ_REDUCED_1, 1, (3, 0), [(4, 0, 1)], id="reduced-one-p2"
    ),
    pytest.param(
        3,
        IntegrationMethod.GGQ_REDUCED_1,
        1,
        (5, 1),
        [(6, 1, 1), (5, 0, 2)],
        id="reduced-one-p3",
    ),
    pytest.param(
        4,
        IntegrationMethod.GGQ_REDUCED_1,
        2,
        (7, 2),
        [(8, 2, 2), (7, 1, 3)],
        id="reduced-one-p4",
    ),
    pytest.param(
        2, IntegrationMethod.GGQ_REDUCED_2, 1, (2, 0), [(3, 0, 2)], id="reduced-two-p2"
    ),
    pytest.param(
        3,
        IntegrationMethod.GGQ_REDUCED_2,
        1,
        (4, 1),
        [(5, 1, 2), (4, 0, 2)],
        id="reduced-two-p3",
    ),
    pytest.param(
        4,
        IntegrationMethod.GGQ_REDUCED_2,
        1,
        (6, 2),
        [(7, 2, 2), (6, 1, 3)],
        id="reduced-two-p4",
    ),
]


@pytest.mark.parametrize(
    ("factory_degree", "method", "start", "passing_rule", "failing_rules"),
    CASES,
)
def test_ggq_rule(
    factory_degree: int,
    method: IntegrationMethod,
    start: int,
    passing_rule: tuple[int, int],
    failing_rules: list[tuple[int, int, int]],
) -> None:
    for num_elements in range(start, 101):
        points = testing.IntegrationPointFactory1D.get_ggq(
            factory_degree, num_elements, method
        )
        _check_ggq_rule(points, *passing_rule, num_elements, True)
        for polynomial_degree, continuity, minimum_elements in failing_rules:
            if num_elements >= minimum_elements:
                _check_ggq_rule(
                    points,
                    polynomial_degree,
                    continuity,
                    num_elements,
                    False,
                )
