"""Tests for public settings normalization and C++ dictionary conversion."""

import json
from collections.abc import Callable
from pathlib import Path

from pyqueso import Model
from pyqueso.scripts.json_io import JsonIO


def test_customized_settings_normalize_to_component_values(
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("settings")
    model = Model(directory / "QuESoSettings_custom_1.json")
    settings = model.settings("main")

    assert settings["component_name"] == "main"
    assert settings["echo_level"] == 2
    assert not settings["write_output_to_file"]
    assert Path(settings["input_filename"]) == directory / "dummy.stl"
    assert Path(settings["output_directory_name"]) == directory / "new_output/main"
    assert settings["background_grid_settings"]["polynomial_order"] == [2, 3, 2]
    assert (
        settings["non_trimmed_quadrature_rule_settings"]["integration_method"]
        == "GGQ_Optimal"
    )
    assert [
        condition["condition_id"] for condition in settings["conditions_settings_list"]
    ] == [1, 2]


def test_defaults_and_json_io_component_conversion(
    copy_test_data: Callable[[str], Path], tmp_path: Path
) -> None:
    directory = copy_test_data("settings")
    model = Model(directory / "QuESoSettings_default.json")
    settings = model.settings("main")
    assert settings["echo_level"] == 1
    assert settings["write_output_to_file"]
    assert Path(settings["output_directory_name"]) == directory / "queso_output/main"

    dictionary = JsonIO.read_settings(settings).dictionary
    assert dictionary.get_string("component_name") == "main"
    assert dictionary.get_int("echo_level") == 1
    assert dictionary["background_grid_settings"].get_int_vector(
        "number_of_elements"
    ) == [1, 1, 1]

    escaped_input = 'C:\\Users\\test\\"quote"\n\x01'
    dictionary.set_value("input_filename", escaped_input)
    filename = tmp_path / "written_settings.json"
    JsonIO.write_settings(dictionary, str(filename))
    written = json.loads(filename.read_text(encoding="utf-8"))
    assert written["component_name"] == "main"
    assert written["echo_level"] == 1
    assert written["input_filename"] == escaped_input
