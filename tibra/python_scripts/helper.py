# Project imports
import TIBRA_PythonApplication as TIBRA_Application

# External imports
import json

def FromPointFromParamToGlobalSpace(point,lower_point, upper_point):
    tmp_point = [0.0, 0.0, 0.0]
    tmp_point[0] = (point[0] * abs(lower_point[0] - upper_point[0])) + lower_point[0]
    tmp_point[1] = (point[1] * abs(lower_point[1] - upper_point[1])) + lower_point[1]
    tmp_point[2] = (point[2] * abs(lower_point[2] - upper_point[2])) + lower_point[2]

    return tmp_point

def FromPointFromGlobalToParamSpace(point,lower_point, upper_point):
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

def ReadParameters(json_filename):
    ''' ReadParameters \n
    Reads Parameters from json file and returns TIBRAApplication.Parameters()
    '''
    with open(json_filename, 'r') as file:
        settings = json.load(file)

    parameters = TIBRA_Application.Parameters()
    items = GetItems(settings)
    for key, value in items:
        parameters.Set(key, value)

    return parameters
