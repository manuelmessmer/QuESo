import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import warnings

class BoundingBox():
    """Base Class to provide interface for the application of boundary conditions.

    Derived class must override 'apply()'.
    """
    def __init__(self, lower_point, upper_point):
        """The constructor."""
        self.lower_point = lower_point
        self.upper_point = upper_point

    def apply(self, model_part):
        raise Exception("BoundingBox :: Base class method is called!")

    def _is_inside(self, point):
        """Simple check if given point is inside the bounding box"""
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
        """Overrides base class method"""
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
        """Overrides base class method"""
        node_id_force = []
        for node in model_part.Nodes:
            point = [node.X0, node.Y0, node.Z0]
            if( self._is_inside(point) ):
                node_id_force.append(node.Id)

        if node_id_force:
            nodal_force_x = self.force[0]/len(node_id_force)
            nodal_force_y = self.force[1]/len(node_id_force)
            nodal_force_z = self.force[2]/len(node_id_force)
            id_counter = id_counter = model_part.NumberOfConditions() + 1
            properties = model_part.GetProperties()[0]
            for i in node_id_force:
                model_part.CreateNewCondition('PointLoadCondition3D1N', id_counter, [i], properties)
                model_part.GetNode(i).SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD_X, nodal_force_x)
                model_part.GetNode(i).SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD_Y, nodal_force_y)
                model_part.GetNode(i).SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD_Z, nodal_force_z)
                id_counter += 1
        else:
            warnings.warn("No nodes are within this boundaing Box", stacklevel=2)

