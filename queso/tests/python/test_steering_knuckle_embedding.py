"""Embedding and component-information regression for the steering knuckle."""

import json
from collections.abc import Callable
from pathlib import Path

import pyqueso
import pytest


def test_steering_knuckle_embedding(
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("steering_knuckle")
    model = pyqueso.Model(directory / "embedding/QuESoSettings.json")
    model.create()

    output = directory / "embedding/output/main/component_info.json"
    output_info = json.loads(output.read_text(encoding="utf-8"))
    represented_volume = sum(
        point.weight
        for element in model.elements("main")
        for point in element.integration_points
    )
    component_info = model.component_info("main")

    quadrature_info = component_info["quadrature_info"]
    assert represented_volume == pytest.approx(
        quadrature_info.get_double("represented_volume"), rel=0.0, abs=5.0e-8
    )

    geometry_info = component_info["embedded_geometry_info"]
    assert geometry_info.get_bool("is_closed")
    assert geometry_info.get_double("volume") == pytest.approx(
        represented_volume, rel=0.0, abs=5.0e-6
    )
    assert output_info["embedded_geometry_info"]["is_closed"]
    assert output_info["embedded_geometry_info"]["volume"] == pytest.approx(
        represented_volume, rel=0.0, abs=5.0e-5
    )

    assert quadrature_info.get_double("percentage_of_geometry_volume") == pytest.approx(
        100.0, rel=0.0, abs=5.0e-6
    )
    assert quadrature_info.get_int("tot_num_points") == 13255
    assert quadrature_info.get_double(
        "num_of_points_per_full_element"
    ) == pytest.approx(23.25, rel=0.0, abs=5.0e-6)
    assert 26 < quadrature_info.get_double("num_of_points_per_trimmed_element") < 27

    output_quadrature = output_info["quadrature_info"]
    assert output_quadrature["percentage_of_geometry_volume"] == pytest.approx(
        100.0, rel=0.0, abs=5.0e-6
    )
    assert output_quadrature["tot_num_points"] == 13255
    assert output_quadrature["num_of_points_per_full_element"] == pytest.approx(
        23.25, rel=0.0, abs=5.0e-6
    )
    assert 26 < output_quadrature["num_of_points_per_trimmed_element"] < 27

    grid_info = component_info["background_grid_info"]
    assert grid_info.get_int("num_active_elements") == 498
    assert grid_info.get_int("num_trimmed_elements") == 486
    assert grid_info.get_int("num_full_elements") == 12
    assert grid_info.get_int("num_inactive_elements") == 5752

    output_grid = output_info["background_grid_info"]
    assert output_grid["num_active_elements"] == 498
    assert output_grid["num_trimmed_elements"] == 486
    assert output_grid["num_full_elements"] == 12
    assert output_grid["num_inactive_elements"] == 5752

    elapsed = component_info["elapsed_time_info"]
    output_elapsed = output_info["elapsed_time_info"]
    assert elapsed.get_double("total") > 0.0
    assert output_elapsed["total"] > 0.0

    volume_time = elapsed["volume_time_info"]
    output_volume_time = output_elapsed["volume_time_info"]
    for name in (
        "total",
        "classification_of_elements",
        "computation_of_intersections",
        "solution_of_moment_fitting_eqs",
        "construction_of_ggq_rules",
    ):
        assert volume_time.get_double(name) > 0.0
        assert output_volume_time[name] > 0.0

    assert elapsed["conditions_time_info"].get_double("total") > 0.0
    assert output_elapsed["conditions_time_info"]["total"] > 0.0
    assert elapsed["write_files_time_info"].get_double("total") > 0.0
    assert output_elapsed["write_files_time_info"]["total"] > 0.0

    conditions_info = component_info.get_list("conditions_infos_list")
    output_conditions = output_info["conditions_infos_list"]
    expected_areas = [332.3775399, 577.9781408, 921.1636352, 1183.543044]
    for index, (info, output_condition, expected_area) in enumerate(
        zip(conditions_info, output_conditions, expected_areas), start=1
    ):
        assert info.get_int("condition_id") == index
        assert output_condition["condition_id"] == index
        assert info.get_double("surf_area") == pytest.approx(
            expected_area, rel=0.0, abs=5.0e-6
        )
        assert output_condition["surf_area"] == pytest.approx(
            expected_area, rel=0.0, abs=5.0e-6
        )
        assert info.get_double("perc_surf_area_in_active_domain") == pytest.approx(
            100.0, rel=0.0, abs=5.0e-6
        )
        assert output_condition["perc_surf_area_in_active_domain"] == pytest.approx(
            100.0, rel=0.0, abs=5.0e-6
        )
