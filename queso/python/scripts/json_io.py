# Project imports
import pyqueso

# External imports
import json

class JsonIO():
    """
    Utility class for reading and writing QuESo settings to and from JSON files.
    """
    @staticmethod
    def write_settings(
            settings: pyqueso.Dictionary, # type: ignore (TODO: add .pyi)
            json_filename: str
        ) -> None:
        """
        Write QuESo settings to a JSON file.

        Args:
            settings (pyqueso.Dictionary): The settings object to write.
            json_filename (str): The path to the output JSON file.
        """
        pyqueso.io.write_dictionary_to_json(settings, json_filename)

    @classmethod
    def read_settings(cls, json_filename: str) -> pyqueso.DictionaryHolder: # type: ignore (TODO: add .pyi)
        """
        Read QuESo settings from a JSON file.

        Args:
            json_filename (str): The path to the JSON file to read.

        Returns:
            pyqueso.DictionaryHolder: QuESo dictionary holder object that contains the parsed settings.
        """
        with open(json_filename, 'r') as file:
            dictionary = json.load(file)

        queso_settings_holder = pyqueso.Dictionary.create("Settings")
        cls._read_dict(dictionary, queso_settings_holder.dictionary)

        return queso_settings_holder

    @classmethod
    def _read_dict(cls,
            dictionary: dict[str, object],
            queso_settings: pyqueso.Dictionary # type: ignore (TODO: add .pyi)
        ) -> None:
        """
        Recursively populate a QuESo settings object from a dictionary.

        Args:
            dictionary (dict[str, object]): Parsed JSON dictionary.
            queso_settings (pyqueso.Dictionary): The settings object to populate.
        """
        for string_key, value in dictionary.items():
            if isinstance(value, dict):
                queso_sub_settings = queso_settings[string_key]
                cls._read_dict( value, queso_sub_settings ) # Got to next level
            elif isinstance(value, list):
                cls._read_list( string_key, value, queso_settings )
            else:
                cls._set_value( string_key, value, queso_settings )

    @classmethod
    def _read_list(cls,
            string_key: str,
            value: list[object],
            queso_settings: pyqueso.Dictionary # type: ignore (TODO: add .pyi)
        ) -> None:
        """
        Read a list-type setting entry.

        Args:
            string_key (str): The setting key.
            value (list): The list of values from JSON.
            queso_settings (pyqueso.Dictionary): The settings object to modify.
        """
        if( isinstance(value[0], dict) ):
            queso_list = queso_settings.get_list(string_key)
            cls._read_conditions_settings_list(value, queso_list)
        else:
            queso_settings.set_value(string_key, value)

    @classmethod
    def _read_conditions_settings_list(cls,
            condition_settings_list: list[dict[str, object]],
            queso_dict_list: pyqueso.DictionaryList # type: ignore (TODO: add .pyi)
        ) -> None:
        """
        Read a list of condition settings from JSON.

        Args:
            condition_settings_list (list[dict[str, object]]): List of condition dictionaries.
            queso_dict_list (pyqueso.DictionaryList): The settings list object to populate.
        """
        for condition_settings in condition_settings_list:
            new_cond_settings_holder = pyqueso.Dictionary.create("ConditionSettings")
            cls._read_dict(condition_settings, new_cond_settings_holder.dictionary)
            queso_dict_list.append(new_cond_settings_holder)


    @classmethod
    def _set_value(cls,
            string_key: str,
            value: object,
            queso_settings: pyqueso.Dictionary # type: ignore (TODO: add .pyi)
        ) -> None:
        """
        Set a single value in the settings, with type conversion for enums.

        Args:
            string_key (str): The setting key.
            value (object): The value from JSON.
            queso_settings (pyqueso.Dictionary): The settings object to modify.
        """
        if( value != "Not Set."):
            if string_key == "integration_method":
                # Convert string to enum
                enum_value = cls._get_enum(value, cls.string_to_enum_integration_method)
                queso_settings.set_value(string_key, enum_value)
            elif string_key == "grid_type":
                # Convert string to enum
                enum_value = cls._get_enum(value, cls.string_to_enum_grid_type)
                queso_settings.set_value(string_key, enum_value)
            else:
                queso_settings.set_value(string_key, value)

    @classmethod
    def _get_enum(cls,
            string_key: str,
            string_to_enum_dict: dict[str, object]
        ) -> object:
        """
        Get enum value from string representation.

        Args:
            string_key (str): The string representation of the enum.
            string_to_enum_dict (dict[str, object]): Mapping from string to enum.

        Returns:
            Enum: The corresponding enum value.

        Raises:
            Exception: If the string_key is not in the mapping.
        """
        if string_key in string_to_enum_dict:
            return string_to_enum_dict[string_key]

        error_msg = (
            f"JsonIO :: Given parameter ({string_key}) not available. "
            f"Possible options: {cls._get_available_keys(string_to_enum_dict)}\n"
        )
        raise Exception(error_msg)

    @classmethod
    def _get_available_keys(cls, string_to_enum_dict: dict[str, object]) -> list[str]:
        """
        Get available keys from a string-to-enum mapping, excluding `_values` suffixes.

        Args:
            string_to_enum_dict (dict[str, object]): Mapping from string to enum.

        Returns:
            list[str]: List of valid keys.
        """
        keys = []
        for key in string_to_enum_dict.keys():
            if not key.endswith("_values"):
                keys.append(key)
        return keys

    string_to_enum_integration_method = {
        "Gauss": pyqueso.IntegrationMethod.GAUSS,
        "Gauss_Reduced1": pyqueso.IntegrationMethod.GAUSS_REDUCED_1,
        "Gauss_Reduced2": pyqueso.IntegrationMethod.GAUSS_REDUCED_2,
        "GGQ_Optimal": pyqueso.IntegrationMethod.GGQ_OPTIMAL,
        "GGQ_Reduced1": pyqueso.IntegrationMethod.GGQ_REDUCED_1,
        "GGQ_Reduced2": pyqueso.IntegrationMethod.GGQ_REDUCED_2,
    }

    string_to_enum_grid_type = {
        "b_spline_grid": pyqueso.GridType.B_SPLINE_GRID,
        "hexahedral_fe_grid": pyqueso.GridType.HEXAHEDRAL_FE_GRID,
    }

