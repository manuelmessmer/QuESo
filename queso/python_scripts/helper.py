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

    return tmp_point

def PointFromParamToGlobalSpace(point, bound_xyz, bound_uvw):
    tmp_point = [0.0, 0.0, 0.0]
    lower_point_xyz = bound_xyz[0]
    upper_point_xyz = bound_xyz[1]
    lower_point_uvw = bound_uvw[0]
    upper_point_uvw = bound_uvw[1]
    tmp_point[0] = (point[0] - lower_point_uvw[0]) * abs(upper_point_xyz[0] - lower_point_xyz[0]) / abs(upper_point_uvw[0]-lower_point_uvw[0]) + lower_point_xyz[0]
    tmp_point[1] = (point[1] - lower_point_uvw[1]) * abs(upper_point_xyz[1] - lower_point_xyz[1]) / abs(upper_point_uvw[1]-lower_point_uvw[1]) + lower_point_xyz[1]
    tmp_point[2] = (point[2] - lower_point_uvw[2]) * abs(upper_point_xyz[2] - lower_point_xyz[2]) / abs(upper_point_uvw[2]-lower_point_uvw[2]) + lower_point_xyz[2]

    return tmp_point

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
    components = GetComponents(settings)

    for key, value in components.items():
        parameters.Set(key, value)

    for condition_id, condition in enumerate(GetConditions(settings)):
        if "dirichlet" in condition:
            value = condition["dirichlet"]
            parameters.AddDirichletBC(condition_id, value["filename"], QuESo_Application.Point(0.0, 0.0, 0.0), value["penalty_factor"])
        elif "neumann" in condition:
            value = condition["neumann"]
            point = QuESo_Application.Point(value["force"])
            parameters.AddNeumannBC(condition_id, value["filename"], point )
        else:
            raise Exception("Helper::ReadParameters :: Condition type does not exits.")

    return parameters
