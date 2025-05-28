from typing import List, Dict, Any

# Project imports
import QuESo_PythonApplication as QuESo

# External imports
import json

class JsonIO():
    """
    Utility class for reading and writing QuESo settings to and from JSON files.
    """
    @staticmethod
    def write_settings(
            settings: QuESo.Dictionary, # type: ignore (TODO: add .pyi)
            json_filename: str
        ) -> None:
        """
        Write QuESo settings to a JSON file.

        Args:
            settings (QuESo.Dictionary): The settings object to write.
            json_filename (str): The path to the output JSON file.
        """
        QuESo.IO.WriteDictionaryToJSON(settings, json_filename) # type: ignore (TODO: add .pyi)

    @classmethod
    def read_settings(cls, json_filename: str) -> QuESo.DictionaryHolder: # type: ignore (TODO: add .pyi)
        """
        Read QuESo settings from a JSON file.

        Args:
            json_filename (str): The path to the JSON file to read.

        Returns:
            QuESo.DictionaryHolder: QuESo dictionary holder object that contains the parsed settings.
        """
        with open(json_filename, 'r') as file:
            dictionary = json.load(file)

        queso_settings_holder = QuESo.Dictionary.Create("Settings") # type: ignore (TODO: add .pyi)
        cls._read_dict(dictionary, queso_settings_holder.Get() )

        return queso_settings_holder

    @classmethod
    def _read_dict(cls,
            dictionary: Dict,
            queso_settings: QuESo.Dictionary # type: ignore (TODO: add .pyi)
        ) -> None:
        """
        Recursively populate a QuESo settings object from a dictionary.

        Args:
            dictionary (Dict): Parsed JSON dictionary.
            queso_settings (QuESo.Dictionary): The settings object to populate.
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
            value: List[Any],
            queso_settings: QuESo.Dictionary # type: ignore (TODO: add .pyi)
        ) -> None:
        """
        Read a list-type setting entry.

        Args:
            string_key (str): The setting key.
            value (list): The list of values from JSON.
            queso_settings (QuESo.Dictionary): The settings object to modify.
        """
        if( isinstance(value[0], dict) ):
            queso_list = queso_settings.GetList(string_key) # Just called to check key, and throw an error if necessary.
            cls._read_conditions_settings_list(value, queso_list)
        else:
            queso_settings.SetValue(string_key, value )

    @classmethod
    def _read_conditions_settings_list(cls,
            condition_settings_list: List[Dict],
            queso_dict_list: QuESo.DictionaryList # type: ignore (TODO: add .pyi)
        ) -> None:
        """
        Read a list of condition settings from JSON.

        Args:
            condition_settings_list (List[Dict]): List of condition dictionaries.
            queso_dict_list (QuESo.DictionaryList): The settings list object to populate.
        """
        for condition_settings in condition_settings_list:
            new_cond_settings_holder = QuESo.Dictionary.Create("ConditionSettings") # type: ignore (TODO: add .pyi)
            cls._read_dict(condition_settings, new_cond_settings_holder.Get())
            queso_dict_list.append(new_cond_settings_holder)


    @classmethod
    def _set_value(cls,
            string_key: str,
            value: Any,
            queso_settings: QuESo.Dictionary # type: ignore (TODO: add .pyi)
        ) -> None:
        """
        Set a single value in the settings, with type conversion for enums.

        Args:
            string_key (str): The setting key.
            value (any): The value from JSON.
            queso_settings (QuESo.Dictionary): The settings object to modify.
        """
        if( value != "Not Set."):
            if string_key == "integration_method":
                # Convert string to enum
                enum_value = cls._get_enum(value, cls.string_to_enum_integration_method)
                queso_settings.SetValue(string_key, enum_value )
            elif string_key == "grid_type":
                # Convert string to enum
                enum_value = cls._get_enum(value, cls.string_to_enum_grid_type)
                queso_settings.SetValue(string_key, enum_value )
            else:
                queso_settings.SetValue(string_key, value )

    @classmethod
    def _get_enum(cls,
            string_key: str,
            string_to_enum_dict: Dict[str, Any]
        ) -> Any:
        """
        Get enum value from string representation.

        Args:
            string_key (str): The string representation of the enum.
            string_to_enum_dict (Dict[str, Any]): Mapping from string to enum.

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
    def _get_available_keys(cls, string_to_enum_dict: Dict[str, Any]) -> Any:
        """
        Get available keys from a string-to-enum mapping, excluding `_values` suffixes.

        Args:
            string_to_enum_dict (Dict[str, Any]): Mapping from string to enum.

        Returns:
            list[Any]: List of valid keys.
        """
        keys = []
        for key in string_to_enum_dict.keys():
            if not key.endswith("_values"):
                keys.append(key)
        return keys

    string_to_enum_integration_method = {
        "Gauss"  : QuESo.IntegrationMethod.Gauss, # type: ignore (TODO: add .pyi)
        "Gauss_Reduced1"  : QuESo.IntegrationMethod.Gauss_Reduced1, # type: ignore (TODO: add .pyi)
        "Gauss_Reduced2"  : QuESo.IntegrationMethod.Gauss_Reduced2, # type: ignore (TODO: add .pyi)
        "GGQ_Optimal"  : QuESo.IntegrationMethod.GGQ_Optimal, # type: ignore (TODO: add .pyi)
        "GGQ_Reduced1"  : QuESo.IntegrationMethod.GGQ_Reduced1, # type: ignore (TODO: add .pyi)
        "GGQ_Reduced2"  : QuESo.IntegrationMethod.GGQ_Reduced2, # type: ignore (TODO: add .pyi)
    }

    string_to_enum_grid_type = {
        "b_spline_grid"  : QuESo.GridType.b_spline_grid, # type: ignore (TODO: add .pyi)
        "hexahedral_fe_grid"  : QuESo.GridType.hexahedral_fe_grid # type: ignore (TODO: add .pyi)
    }

