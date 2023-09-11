# Project imports
import QuESo_PythonApplication as QuESo_Application

# External imports
import json


def PointFromGlobalToParamSpace(point, bound_xyz, bound_uvw):
    tmp_point = [0.0, 0.0, 0.0]
    lower_point_xyz = bound_xyz[0]
    upper_point_xyz = bound_xyz[1]
    lower_point_uvw = bound_uvw[0]
    upper_point_uvw = bound_uvw[1]
    tmp_point[0] = (point[0] - lower_point_xyz[0]) * abs(upper_point_uvw[0] - lower_point_uvw[0]) / abs(upper_point_xyz[0]-lower_point_xyz[0]) + lower_point_uvw[0]
    tmp_point[1] = (point[1] - lower_point_xyz[1]) * abs(upper_point_uvw[1] - lower_point_uvw[1]) / abs(upper_point_xyz[1]-lower_point_xyz[1]) + lower_point_uvw[1]
    tmp_point[2] = (point[2] - lower_point_xyz[2]) * abs(upper_point_uvw[2] - lower_point_uvw[2]) / abs(upper_point_xyz[2]-lower_point_xyz[2]) + lower_point_uvw[2]

    return QuESo_Application.Point(tmp_point[0], tmp_point[1], tmp_point[2])

def PointFromParamToGlobalSpace(point, bound_xyz, bound_uvw):
    tmp_point = [0.0, 0.0, 0.0]
    lower_point_xyz = bound_xyz[0]
    upper_point_xyz = bound_xyz[1]
    lower_point_uvw = bound_uvw[0]
    upper_point_uvw = bound_uvw[1]
    tmp_point[0] = (point[0] - lower_point_uvw[0]) * abs(upper_point_xyz[0] - lower_point_xyz[0]) / abs(upper_point_uvw[0]-lower_point_uvw[0]) + lower_point_xyz[0]
    tmp_point[1] = (point[1] - lower_point_uvw[1]) * abs(upper_point_xyz[1] - lower_point_xyz[1]) / abs(upper_point_uvw[1]-lower_point_uvw[1]) + lower_point_xyz[1]
    tmp_point[2] = (point[2] - lower_point_uvw[2]) * abs(upper_point_xyz[2] - lower_point_xyz[2]) / abs(upper_point_uvw[2]-lower_point_uvw[2]) + lower_point_xyz[2]

    return QuESo_Application.Point(tmp_point[0], tmp_point[1], tmp_point[2])

def GetItems(dictionary):
    for k, v in dictionary.items():
        if isinstance(v, dict):
            yield from GetItems(v)
        else:
            yield (k, v)


def GetComponents(dictionary):
    ''' Extract components from dictionary. If 'conditions' is not available, return empty dictionary.
    '''
    components = {}
    # Get required setttings
    if "general_settings" in dictionary:
        components.update(GetItems(dictionary["general_settings"]))
    else:
        raise Exception("Helper::GetComponents :: Json does not contain 'general_settings'.")

    if "mesh_settings" in dictionary:
        components.update(GetItems(dictionary["mesh_settings"]))
    else:
        raise Exception("Helper::GetComponents :: Json does not contain 'mesh_settings'.")

    # Get non-required setttings
    if "trimmed_quadrature_rule_settings" in dictionary:
        components.update(GetItems(dictionary["trimmed_quadrature_rule_settings"]))

    if "non_trimmed_quadrature_rule_settings" in dictionary:
        components.update(GetItems(dictionary["non_trimmed_quadrature_rule_settings"]))

    return components

def GetConditions(dictionary):
    ''' Extract conditions from dictionary. If 'conditions' is not available, return empty dictionary.
    '''
    if "conditions" in dictionary:
        return dictionary["conditions"]
    else:
        return {}

def ReadParameters(json_filename):
    ''' ReadParameters \n
    Reads Parameters from json file and returns QuESoApplication.Parameters()
    parameters = QuESo_Application.Parameters()
    '''
    with open(json_filename, 'r') as file:
        settings = json.load(file)

    parameters = QuESo_Application.Parameters()

    # Read general settings
    general_setting = QuESo_Application.GlobalParameters()
    for key, value in GetComponents(settings).items():
        general_setting.Set(key, value)

    parameters.AddGlobalParameters(general_setting)

    # Read general settings
    for condition in GetConditions(settings):
        type = list(condition.keys())[0]
        condition_settings = QuESo_Application.ConditionParameters(str(type))
        for key, value in condition[type].items():
            condition_settings.Set(key, value)

        parameters.AddConditionParameters(condition_settings)

    return parameters
