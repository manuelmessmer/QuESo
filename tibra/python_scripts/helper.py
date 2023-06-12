# Project imports
import TIBRA_PythonApplication as TIBRA_Application

# External imports
import json

def PointFromParamToGlobalSpace(point,lower_point, upper_point):
    tmp_point = [0.0, 0.0, 0.0]
    tmp_point[0] = (point[0] * abs(lower_point[0] - upper_point[0])) + lower_point[0]
    tmp_point[1] = (point[1] * abs(lower_point[1] - upper_point[1])) + lower_point[1]
    tmp_point[2] = (point[2] * abs(lower_point[2] - upper_point[2])) + lower_point[2]

    return tmp_point

def PointFromGlobalToParamSpace(point,lower_point, upper_point):
    tmp_point = [0.0, 0.0, 0.0]
    tmp_point[0] = (point[0] - lower_point[0])/ abs(lower_point[0] - upper_point[0])
    tmp_point[1] = (point[1] - lower_point[1])/ abs(lower_point[1] - upper_point[1])
    tmp_point[2] = (point[2] - lower_point[2])/ abs(lower_point[2] - upper_point[2])

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
    Reads Parameters from json file and returns TIBRAApplication.Parameters()
    parameters = TIBRA_Application.Parameters()
    '''
    with open(json_filename, 'r') as file:
        settings = json.load(file)

    parameters = TIBRA_Application.Parameters()
    components = GetComponents(settings)

    for key, value in components.items():
        parameters.Set(key, value)

    for condition_id, condition in enumerate(GetConditions(settings)):
        if "dirichlet" in condition:
            value = condition["dirichlet"]
            #print(value)
            parameters.AddDirichletBC(condition_id, value["filename"], TIBRA_Application.Point(0.0, 0.0, 0.0), value["penalty_factor"])
        elif "neumann" in condition:
            value = condition["neumann"]
            #print(value)
            point = TIBRA_Application.Point(value["force"])
            parameters.AddNeumannBC(condition_id, value["filename"], point )
        else:
            raise Exception("Helper::ReadParameters :: Condition type does not exits.")

    return parameters
