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
        if k != "conditions":
            if isinstance(v, dict):
                yield from GetItems(v)
            else:
                yield (k, v)

def GetConditions(dictionary):
    for k, v in dictionary.items():
        if k == "conditions":
            yield from GetConditions(v)
        if k == "neumann" or k == "dirichlet":
            yield (k, v)

# def GetDirichletConditions(dictionary):
#     for k, v in dictionary.items():
#         if k == "Conditions":
#             yield from GetItems(v)
#         if k == "Dirichlet":
#             yield (k, v)


def ReadParameters(json_filename):
    ''' ReadParameters \n
    Reads Parameters from json file and returns TIBRAApplication.Parameters()
    '''
    with open(json_filename, 'r') as file:
        settings = json.load(file)

    parameters = TIBRA_Application.Parameters()
    items = GetItems(settings)
    for key, value in items:
        print("Key: ", key)
        print("Value: ", value)
        parameters.Set(key, value)

    for key, value in GetConditions(settings):
        if key == "dirichlet":
            point = TIBRA_Application.Point(value["displacement"])
            parameters.AddCondition("Dirichlet", value["filename"], point, value["penalty_factor"])
        if key == "neumann":
            point = TIBRA_Application.Point(value["force"])
            parameters.AddCondition("Dirichlet", value["filename"], point, 0.0)

    return parameters
