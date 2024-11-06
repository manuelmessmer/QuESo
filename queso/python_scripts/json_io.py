# Project imports
import QuESo_PythonApplication as QuESo_Application

# External imports
import json

class JsonIO():
    @classmethod
    def WriteSettings(cls, settings, json_filename):
        """ Write settings to json file.

        @param settings
        @param json_filename
        """
        QuESo_Application.IO.WriteSettingsToJSON(settings, json_filename)

    @classmethod
    def ReadSettings(cls, json_filename):
        ''' Reads settings from json file and returns QuESoApplication.Settings()

        @param json_filename

        @return QuESo_Application.Settings
        '''
        with open(json_filename, 'r') as file:
            dictionary = json.load(file)

        queso_settings = QuESo_Application.Settings()
        cls._ReadDict(dictionary, queso_settings )

        return queso_settings

    @classmethod
    def _ReadDict(cls, dictionary, queso_settings ):
        for string_key, value in dictionary.items():
            string_key_for_values = string_key + "_values"
            if isinstance(value, dict):
                queso_sub_settings = queso_settings[string_key]
                cls._ReadDict( value, queso_sub_settings ) # Got to next level
            elif isinstance(value, list):
                cls._ReadList( string_key, value, queso_settings )
            else:
                cls._SetValue( string_key, value, queso_settings )

    @classmethod
    def _ReadList(cls, string_key, value, queso_settings):
        if( isinstance(value[0], dict) ):
            queso_settings.GetList(string_key) # Just called to check key, and throw an error if necessary.
            cls._ReadConditionsSettingsList(value, queso_settings)
        else:
            queso_settings.SetValue(string_key, value )

    @classmethod
    def _ReadConditionsSettingsList(cls, condition_settings_list, queso_settings):
        for condition_settings in condition_settings_list:
            new_cond_settings = queso_settings.CreateNewConditionSettings()
            cls._ReadDict(condition_settings, new_cond_settings)

    @classmethod
    def _SetValue(cls, string_key, value, queso_settings):
        if( value != "Not Set."):
            if string_key == "integration_method":
                # Convert string to enum
                enum_value = cls._GetEnum(value, cls.string_to_enum_integration_method)
                queso_settings.SetValue(string_key, enum_value )
            elif string_key == "grid_type":
                # Convert string to enum
                enum_value = cls._GetEnum(value, cls.string_to_enum_grid_type)
                queso_settings.SetValue(string_key, enum_value )
            else:
                queso_settings.SetValue(string_key, value )

    @classmethod
    def _GetEnum(cls, string_key, string_to_enum_dict):
        if string_key in string_to_enum_dict:
            return string_to_enum_dict[string_key]

        error_msg = "JsonIO :: Given parameter (" + string_key + ") not available. Possible options: " + str(cls._GetAvailableKeys(string_to_enum_dict)) + '\n'
        raise Exception(error_msg)

    @classmethod
    def _GetAvailableKeys(cls, string_to_enum_dict):
        keys = []
        for key in string_to_enum_dict.keys():
            if not key.endswith("_values"):
                keys.append(key)
        return keys

    string_to_enum_integration_method = {
        "Gauss"  : QuESo_Application.IntegrationMethod.Gauss,
        "Gauss_Reduced1"  : QuESo_Application.IntegrationMethod.Gauss_Reduced1,
        "Gauss_Reduced2"  : QuESo_Application.IntegrationMethod.Gauss_Reduced2,
        "GGQ_Optimal"  : QuESo_Application.IntegrationMethod.GGQ_Optimal,
        "GGQ_Reduced1"  : QuESo_Application.IntegrationMethod.GGQ_Reduced1,
        "GGQ_Reduced2"  : QuESo_Application.IntegrationMethod.GGQ_Reduced2,
    }

    string_to_enum_grid_type = {
        "b_spline_grid"  : QuESo_Application.GridType.b_spline_grid,
        "hexahedral_fe_grid"  : QuESo_Application.GridType.hexahedral_fe_grid
    }

