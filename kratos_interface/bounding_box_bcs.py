from typing import Tuple
import warnings

# Kratos imports
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Type definition
Point3D = Tuple[float, float, float]

class BoundingBoxBC():
    """
    Base class to provide interface for the application of boundary conditions via bounding boxes.

    Derived classes must override the 'apply()' method to implement specific boundary condition logic.
    """
    def __init__(self,
            lower_point: Point3D,
            upper_point: Point3D
        ) -> None:
        """
        Constructor for the BoundingBoxBC class.

        Args:
            lower_point (Point3D): The lower corner of the bounding box (e.g., [x_min, y_min, z_min]).
            upper_point (Point3D): The upper corner of the bounding box (e.g., [x_max, y_max, z_max]).
        """
        self.lower_point = lower_point
        self.upper_point = upper_point

    def apply(self, model_part: KM.ModelPart) -> None:
        """
        Applies the boundary condition to the given model part.

        Args:
            model_part (KM.ModelPart): The model part to apply the condition to.

        Raises:
            Exception: Raises an exception because this method must be overridden in derived classes.
        """
        raise Exception("BoundingBoxBC :: Base class method is called!")

    def _is_inside(self, point: Point3D) -> bool:
        """
        Checks if a given point lies within the defined bounding box.

        Args:
            point (Point3D): A 3D point specified as a tuple [x, y, z].

        Returns:
            bool: True if the point lies within the bounding box; False otherwise.
        """
        return all(self.lower_point[i] < point[i] < self.upper_point[i] for i in range(3))

    @staticmethod
    def is_weak_condition() -> bool:
        """
        Indicates whether the boundary condition is weak.

        Returns:
            bool: Always False.
        """
        return False

class DirichletCondition(BoundingBoxBC):
    """
    Dirichlet boundary condition applied to a given bounding box.

    This class overrides the apply method from the BoundingBoxBC class to
    fix the displacement of nodes inside the bounding box, based on the
    provided displacement vector.
    """
    def __init__(self,
            lower_point: Point3D,
            upper_point: Point3D,
            disp: Tuple[bool, bool, bool]
        ) -> None:
        """
        Initializes the DirichletCondition with bounding box corners and displacement vector.

        Args:
            lower_point (Point3D): The lower corner of the bounding box (e.g., [x_min, y_min, z_min]).
            upper_point (Point3D): The upper corner of the bounding box (e.g., [x_max, y_max, z_max]).
            disp (Tuple[bool, bool, bool]): A vector indicating the displacement conditions for each axis.
                (e.g., [True, False, True] would fix displacement in the X and Z directions).
        """
        super(DirichletCondition, self).__init__(lower_point, upper_point)
        self.disp = disp

    def apply(self, model_part: KM.ModelPart) -> None:
        """
        Applies the Dirichlet boundary condition to the given model part by fixing displacements
        of nodes that lie within the bounding box.

        Args:
            model_part (KM.ModelPart): Kratos model part.
        """
        for node in model_part.Nodes:
            if( self._is_inside([node.X0, node.Y0, node.Z0]) ):
                if self.disp[0]:
                    node.Fix(KM.DISPLACEMENT_X)
                if self.disp[1]:
                    node.Fix(KM.DISPLACEMENT_Y)
                if self.disp[2]:
                    node.Fix(KM.DISPLACEMENT_Z)

class NeumannCondition(BoundingBoxBC):
    """
    Neumann boundary condition applied to a given bounding box.

    This class overrides the `apply` method from the `BoundingBoxBC` class to apply a force
    to the nodes inside the bounding box. The applied force is distributed equally among
    the nodes within the bounding box.
    """
    def __init__(self,
            lower_point: Point3D,
            upper_point: Point3D,
            force: Tuple[float, float, float]
        ) -> None:
        """
        Initializes the NeumannCondition with bounding box corners and a force vector.

        Args:
            lower_point (Point3D): The lower corner of the bounding box (e.g., [x_min, y_min, z_min]).
            upper_point (Point3D): The upper corner of the bounding box (e.g., [x_max, y_max, z_max]).
            force (Tuple[float, float, float]): The total force to be applied in the X, Y, and Z directions.
                (e.g., [fx, fy, fz]).
        """
        super(NeumannCondition, self).__init__(lower_point, upper_point)
        self.force = force

    def apply(self, model_part: KM.ModelPart) -> None:
        """
        Applies the Neumann boundary condition to the given model part by distributing
        the force across the nodes inside the bounding box.

        Args:
            model_part (KM.ModelPart): The model or part of the model to apply the boundary condition to.

        Warnings:
            If no nodes are inside the bounding box, a warning is raised.
        """

        node_ids_with_force = [node.Id for node in model_part.Nodes if self._is_inside([node.X0, node.Y0, node.Z0])]

        if node_ids_with_force:
            num_nodes = len(node_ids_with_force)
            nodal_force_x, nodal_force_y, nodal_force_z = [f / num_nodes for f in self.force]

            id_counter = id_counter = model_part.NumberOfConditions() + 1
            properties = model_part.GetProperties()[1]
            for i, node_id in enumerate(node_ids_with_force):
                model_part.CreateNewCondition('PointLoadCondition3D1N', id_counter+i, [node_id], properties)
                model_part.GetNode(node_id).SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD_X, nodal_force_x)
                model_part.GetNode(node_id).SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD_Y, nodal_force_y)
                model_part.GetNode(node_id).SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD_Z, nodal_force_z)
        else:
            warnings.warn("There are no nodes within this bounding box.")

