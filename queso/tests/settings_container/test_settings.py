"""Tests for public settings normalization and C++ dictionary conversion."""

import json
import tempfile
import unittest
from pathlib import Path

import pyqueso
from pyqueso.scripts.json_io import JsonIO

DIRECTORY = Path(__file__).parent


class TestSettingsContainer(unittest.TestCase):
    """Exercise public settings fixtures through the multi-component API."""

    def test_customized_settings_normalize_to_component_values(self) -> None:
        model = pyqueso.Model(DIRECTORY / "QuESoSettings_custom_1.json")
        settings = model.settings("main")

        self.assertEqual(settings["component_name"], "main")
        self.assertEqual(settings["echo_level"], 2)
        self.assertFalse(settings["write_output_to_file"])
        self.assertEqual(Path(settings["input_filename"]), DIRECTORY / "dummy.stl")
        self.assertEqual(
            Path(settings["output_directory_name"]), DIRECTORY / "new_output/main"
        )
        self.assertEqual(
            settings["background_grid_settings"]["polynomial_order"], [2, 3, 2]
        )
        self.assertEqual(
            settings["non_trimmed_quadrature_rule_settings"]["integration_method"],
            "GGQ_Optimal",
        )
        self.assertEqual(
            [
                condition["condition_id"]
                for condition in settings["conditions_settings_list"]
            ],
            [1, 2],
        )

    def test_defaults_and_json_io_component_conversion(self) -> None:
        model = pyqueso.Model(DIRECTORY / "QuESoSettings_default.json")
        settings = model.settings("main")
        self.assertEqual(settings["echo_level"], 1)
        self.assertTrue(settings["write_output_to_file"])
        self.assertEqual(
            Path(settings["output_directory_name"]), DIRECTORY / "queso_output/main"
        )

        holder = JsonIO.read_settings(settings)
        dictionary = holder.dictionary
        self.assertEqual(dictionary.get_string("component_name"), "main")
        self.assertEqual(dictionary.get_int("echo_level"), 1)
        self.assertEqual(
            dictionary["background_grid_settings"].get_int_vector("number_of_elements"),
            [1, 1, 1],
        )

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "settings.json"
            JsonIO.write_settings(dictionary, str(filename))
            written = json.loads(filename.read_text(encoding="utf-8"))
        self.assertEqual(written["component_name"], "main")
        self.assertEqual(written["echo_level"], 1)


if __name__ == "__main__":
    unittest.main()
