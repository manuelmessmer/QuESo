import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

class BoundingBox():
    def __init__(self, lower_point, upper_point):
        self.lower_point = lower_point
        self.upper_point = upper_point

    def apply(self, model_part):
        raise Exception("BoundingBox :: Function of base class is called!")

    def _is_inside(self, point):
        x_value =  point[0] < self.upper_point[0] and point[0] > self.lower_point[0]
        y_value =  point[1] < self.upper_point[1] and point[1] > self.lower_point[1]
        z_value =  point[2] < self.upper_point[2] and point[2] > self.lower_point[2]

        if x_value and y_value and z_value:
            return True
        else:
            return False

class DirichletCondition(BoundingBox):
    def __init__(self, lower_point, upper_point, disp):
        super(DirichletCondition, self).__init__(lower_point, upper_point)
        self.disp = disp

    def apply(self, model_part):
        for node in model_part.Nodes:
            point = [node.X0, node.Y0, node.Z0]
            if( self._is_inside(point) ):
                if self.disp[0] == 1:
                    node.Fix(KM.DISPLACEMENT_X)
                if self.disp[1] == 1:
                    node.Fix(KM.DISPLACEMENT_Y)
                if self.disp[2] == 1:
                    node.Fix(KM.DISPLACEMENT_Z)

class NeumannCondition(BoundingBox):
    def __init__(self, lower_point, upper_point, force):
        super(NeumannCondition, self).__init__(lower_point, upper_point)
        self.force = force

    def apply(self, model_part):
        node_id_force = []
        for node in model_part.Nodes:
            point = [node.X0, node.Y0, node.Z0]
            if( self._is_inside(point) ):
                node_id_force.append(node.Id)

        nodal_force = self.force/len(node_id_force)
        id_counter = id_counter = model_part.NumberOfConditions() + 1
        properties = model_part.GetProperties()[0]
        for i in node_id_force:
            model_part.CreateNewCondition('PointLoadCondition3D1N', id_counter, [i], properties)
            model_part.GetNode(i).SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD_Y, nodal_force)
            id_counter = id_counter + 1
