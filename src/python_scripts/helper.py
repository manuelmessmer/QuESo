
def FromParamToGlobalSpace(point,lower_point, upper_point):
    tmp_point = [0.0, 0.0, 0.0]
    tmp_point[0] = (point[0] * abs(lower_point[0] - upper_point[0])) + lower_point[0]
    tmp_point[1] = (point[1] * abs(lower_point[1] - upper_point[1])) + lower_point[1]
    tmp_point[2] = (point[2] * abs(lower_point[2] - upper_point[2])) + lower_point[2]

    return tmp_point

def FromGlobalToParamSpace(point,lower_point, upper_point):
    tmp_point = [0.0, 0.0, 0.0]
    tmp_point[0] = (point[0] - lower_point[0])/ abs(lower_point[0] - upper_point[0])
    tmp_point[1] = (point[1] - lower_point[1])/ abs(lower_point[1] - upper_point[1])
    tmp_point[2] = (point[2] - lower_point[2])/ abs(lower_point[2] - upper_point[2])

    return tmp_point
