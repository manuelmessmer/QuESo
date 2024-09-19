# Project imports
import QuESo_PythonApplication as QuESo_Application

# External imports
import json

class JsonImport():
    @classmethod
    def ReadSettings(cls, json_filename):
        ''' Reads settings from json file and returns QuESoApplication.Settings()

        @param json_filename

        @return QuESo_Application.Settings
        '''
        with open(json_filename, 'r') as file:
            dictionary = json.load(file)

        queso_settings = QuESo_Application.Settings()
        cls._ReadDict(dictionary, queso_settings, cls.string_to_enum_dict)

        return queso_settings

    @classmethod
    def _ReadDict(cls, dictionary, queso_settings, string_to_enum_dict):
        for string_key, value in dictionary.items():
            string_key_for_values = string_key + "_values"
            if isinstance(value, dict):
                enum_key = cls._GetEnum(string_key, string_to_enum_dict)
                queso_sub_settings = queso_settings[enum_key]
                string_to_emum_sub_dict = string_to_enum_dict[string_key_for_values]
                cls._ReadDict( value, queso_sub_settings, string_to_emum_sub_dict ) # Got to next level
            elif isinstance(value, list):
                enum_key = cls._GetEnum(string_key, string_to_enum_dict)
                cls._ReadList( enum_key, value, queso_settings )
            else:
                enum_key = cls._GetEnum(string_key, string_to_enum_dict)
                cls._SetValue( enum_key, value, queso_settings )

    @classmethod
    def _ReadList(cls, enum_key, value, queso_settings):
        if( enum_key is  QuESo_Application.MainSettings.conditions_settings_list ):
            cls._ReadConditionsSettingsList(value, queso_settings)
        else:
            queso_settings.SetValue(enum_key, value )

    @classmethod
    def _ReadConditionsSettingsList(cls, condition_settings_list, queso_settings):
        for condition_settings in condition_settings_list:
            new_cond_settings = queso_settings.CreateNewConditionSettings()
            string_to_enum_sub_dict = cls.string_to_enum_dict["conditions_settings_list_values"]
            for string_key, value in condition_settings.items():
                enum_key = cls._GetEnum(string_key, string_to_enum_sub_dict)
                cls._SetValue(enum_key, value, new_cond_settings )

    @classmethod
    def _SetValue(cls, enum_key, value, queso_settings):
        if enum_key is QuESo_Application.NonTrimmedQuadratureRuleSettings.integration_method:
            # In this case, also value is of type enum
            enum_value = cls._GetEnum(value, cls.string_to_enum_integration_method)
            queso_settings.SetValue(enum_key, enum_value )
        elif enum_key is QuESo_Application.BackgroundGridSettings.grid_type:
            # In this case, also value is of type enum
            enum_value = cls._GetEnum(value, cls.string_to_enum_grid_type)
            queso_settings.SetValue(enum_key, enum_value )
        else:
            queso_settings.SetValue(enum_key, value )

    @classmethod
    def _GetEnum(cls, string_key, string_to_enum_dict):
        if string_key in string_to_enum_dict:
            return string_to_enum_dict[string_key]

        error_msg = "JsonImport :: Given parameter (" + string_key + ") not available. Possible options: " + str(cls._GetAvailableKeys(string_to_enum_dict)) + '\n'
        raise Exception(error_msg)

    @classmethod
    def _GetAvailableKeys(cls, string_to_enum_dict):
        keys = []
        for key in string_to_enum_dict.keys():
            if not key.endswith("_values"):
                keys.append(key)
        return keys

    string_to_enum_dict = {
        "general_settings"  : QuESo_Application.MainSettings.general_settings,
        "general_settings_values"  : {
            "input_filename"  : QuESo_Application.GeneralSettings.input_filename,
            "output_directory_name"  : QuESo_Application.GeneralSettings.output_directory_name,
            "echo_level"  : QuESo_Application.GeneralSettings.echo_level,
            "write_output_to_file" : QuESo_Application.GeneralSettings.write_output_to_file
        },
        "background_grid_settings"  : QuESo_Application.MainSettings.background_grid_settings,
        "background_grid_settings_values"  : {
            "grid_type"        : QuESo_Application.BackgroundGridSettings.grid_type,
            "lower_bound_xyz"  : QuESo_Application.BackgroundGridSettings.lower_bound_xyz,
            "upper_bound_xyz"  : QuESo_Application.BackgroundGridSettings.upper_bound_xyz,
            "lower_bound_uvw"  : QuESo_Application.BackgroundGridSettings.lower_bound_uvw,
            "upper_bound_uvw"  : QuESo_Application.BackgroundGridSettings.upper_bound_uvw,
            "polynomial_order"  : QuESo_Application.BackgroundGridSettings.polynomial_order,
            "number_of_elements"  : QuESo_Application.BackgroundGridSettings.number_of_elements,
        },
        "trimmed_quadrature_rule_settings"  : QuESo_Application.MainSettings.trimmed_quadrature_rule_settings,
        "trimmed_quadrature_rule_settings_values"  : {
            "moment_fitting_residual"  : QuESo_Application.TrimmedQuadratureRuleSettings.moment_fitting_residual,
            "min_element_volume_ratio"  : QuESo_Application.TrimmedQuadratureRuleSettings.min_element_volume_ratio,
            "min_num_boundary_triangles"  : QuESo_Application.TrimmedQuadratureRuleSettings.min_num_boundary_triangles,
            "neglect_elements_if_stl_is_flawed" : QuESo_Application.TrimmedQuadratureRuleSettings.neglect_elements_if_stl_is_flawed
        },
        "non_trimmed_quadrature_rule_settings"  : QuESo_Application.MainSettings.non_trimmed_quadrature_rule_settings,
        "non_trimmed_quadrature_rule_settings_values"  : {
            "integration_method"  : QuESo_Application.NonTrimmedQuadratureRuleSettings.integration_method
        },
        "conditions_settings_list" : QuESo_Application.MainSettings.conditions_settings_list,
        "conditions_settings_list_values" : {
            "condition_id"  :   QuESo_Application.ConditionSettings.condition_id,
            "condition_type"  :   QuESo_Application.ConditionSettings.condition_type,
            "input_filename"  :   QuESo_Application.ConditionSettings.input_filename,
            "modulus"  :   QuESo_Application.ConditionSettings.modulus,
            "direction"  :   QuESo_Application.ConditionSettings.direction,
            "value"  :   QuESo_Application.ConditionSettings.value,
            "penalty_factor"  :   QuESo_Application.ConditionSettings.penalty_factor
        }
    }

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

