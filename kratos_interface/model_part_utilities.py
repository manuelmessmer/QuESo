from typing import Literal, Tuple
# Import Kratos
import KratosMultiphysics as KM
# Import QuESo modules
import QuESoPythonModule as QuESo
from QuESoPythonModule.kratos_interface.weak_bcs import PenaltySupport
from QuESoPythonModule.kratos_interface.weak_bcs import LagrangeSupport
from QuESoPythonModule.kratos_interface.weak_bcs import SurfaceLoad
from QuESoPythonModule.kratos_interface.weak_bcs import PressureLoad

# Type definition
Point3D = Tuple[float, float, float]

class ModelPartUtilities:
    """
    Utility class for managing Kratos ModelPart operations related to QuESo triangle meshes, elements, and conditions.

    This class provides various static methods to interact with Kratos model parts, including reading from and writing
    to STL files, transferring triangle meshes between QuESo and Kratos, and adding elements and conditions to Kratos
    model parts.
    """

    @staticmethod
    def _write_model_part_to_stl(
            KratosModelPart: KM.ModelPart,
            Filename: str
        ) -> None:
        """
        Deprecated method. Writes a Kratos ModelPart to an STL file.

        Args:
            KratosModelPart (KM.ModelPart): The Kratos ModelPart to be written to the STL file.
            Filename (str): The path to the STL file to be written.

        Deprecated:
            This method is deprecated.
        """
        io_settings = KM.Parameters(
            """{
                "open_mode" : "write",
                "new_entity_type" : "geometry"
            }""")

        stl_io = KM.StlIO(Filename, io_settings) # type: ignore
        stl_io.WriteModelPart(KratosModelPart)

    @staticmethod
    def read_model_part_from_triangle_mesh(
            KratosModelPart: KM.ModelPart,
            TriangleMesh: QuESo.TriangleMesh # type: ignore (TODO: add .pyi)
        ) -> None:
        """
        Reads a Kratos ModelPart from a QuESo triangle mesh.

        Args:
            KratosModelPart (KM.ModelPart): The Kratos ModelPart to populate with data from the triangle mesh.
            TriangleMesh (QuESo.TriangleMesh): The QuESo triangle mesh from which the model part is read.
        """
        vertices = TriangleMesh.GetVertices()
        for v_id, vertex in enumerate(vertices):
            KratosModelPart.CreateNewNode(v_id+1, vertex[0], vertex[1], vertex[2])

        triangles = TriangleMesh.GetTriangles()
        for t_id, triangle in enumerate(triangles):
            KratosModelPart.CreateNewElement("ShellThinElement3D3N", t_id+1, [triangle[0]+1, triangle[1]+1, triangle[2]+1], KratosModelPart.GetProperties()[1])

    @staticmethod
    def read_triangle_mesh_grom_model_part(
            TriangleMesh: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            KratosModelPart: KM.ModelPart,
            type: Literal["Elements", "Conditions"] ="Elements"
        ) -> None:
        """
        Reads a QuESo triangle mesh from a Kratos ModelPart.

        Args:
            TriangleMesh (QuESo.TriangleMesh): The QuESo triangle mesh to populate.
            KratosModelPart (KM.ModelPart): The Kratos ModelPart containing the elements or conditions to be read.
            type (Literal["Elements", "Conditions"]): Type of entities to read from the model part. Options are "Elements" or "Conditions". Defaults to "Elements".

        Raises:
            Exception: If the type is not one of the allowed options or if the geometry does not consist of triangles.
        """
        id_map = {node.Id: queso_id for queso_id, node in enumerate(KratosModelPart.Nodes)}
        for node in KratosModelPart.Nodes: # type: ignore
            TriangleMesh.AddVertex([node.X, node.Y, node.Z])

        if( type == "Elements"):
            entity_list = KratosModelPart.Elements
        elif(type == "Conditions"):
            entity_list = KratosModelPart.Conditions
        else:
            message = (
                f"Given condition type '{type_name}' is not available. "
                f"Available options are: 'Elements' and 'Conditions'.")
            raise Exception(message)

        for entity in entity_list: # type: ignore
            geometry = entity.GetGeometry()

            # Ensure the geometry is a triangle
            if( len(geometry) != 3 ):
                raise Exception("ModelPartUtilities :: read_triangle_mesh_grom_model_part :: Queso only allows triangles.")

            # Get the node ids for the triangle
            node_ids = [id_map[node.Id] for node in geometry]
            TriangleMesh.AddTriangle(node_ids)

            # Compute the normal for the triangle
            normal = QuESo.Triangle.NormalStatic( # type: ignore (TODO: add .pyi)
                [node.X, node.Y, node.Z] for node in geometry
            )
            TriangleMesh.AddNormal(normal)

    @staticmethod
    def _create_bc_from_settings(condition_settings, segment, nurbs_volume_name: str = "NurbsVolume"):
        type_name = condition_settings["condition_type"]
        if type_name == "PenaltySupportCondition":
            return PenaltySupport(
                segment.GetTriangleMesh(),
                None,
                None,
                condition_settings["value"],
                condition_settings["penalty_factor"],
                nurbs_volume_name,
            )
        if type_name == "LagrangeSupportCondition":
            return LagrangeSupport(
                segment.GetTriangleMesh(),
                None,
                None,
                condition_settings["value"],
                nurbs_volume_name,
            )
        if type_name == "SurfaceLoadCondition":
            return SurfaceLoad(
                segment.GetTriangleMesh(),
                None,
                None,
                condition_settings["modulus"],
                condition_settings["direction"],
                nurbs_volume_name,
            )
        if type_name == "PressureLoadCondition":
            return PressureLoad(
                segment.GetTriangleMesh(),
                None,
                None,
                condition_settings["modulus"],
                nurbs_volume_name,
            )
        available_conditions = (
            "PenaltySupportCondition",
            "LagrangeSupportCondition",
            "SurfaceLoadCondition",
            "PressureLoadCondition",
        )
        options = ", ".join(f'"{cond}"' for cond in available_conditions)
        raise Exception(f"Given condition type '{type_name}' is not available. Available options are: {options}.")

    @staticmethod
    def add_elements_to_model_part(
            KratosNurbsVolumeModelPart: KM.ModelPart,
            Elements: QuESo.ElementVector, # type: ignore (TODO: add .pyi)
            properties_id: int,
            nurbs_volume_name: str = "NurbsVolume",
            first_element_id: int | None = None,
            element_name: str = "UpdatedLagrangianElement3D8N"
        ) -> int:
        """
        Adds the QuESo elements to the Kratos NurbsVolume ModelPart.

        Args:
            KratosNurbsVolumeModelPart (KM.ModelPart): The Kratos model part to which the elements will be added.
            Elements (QuESo.ElementVector): List of QuESo elements to add to the model part.
        """
        nurbs_volume = KratosNurbsVolumeModelPart.GetGeometry(nurbs_volume_name)
        volume_properties = KratosNurbsVolumeModelPart.GetProperties()[properties_id]

        if first_element_id is None:
            el_count = KratosNurbsVolumeModelPart.NumberOfElements()
        else:
            el_count = int(first_element_id) - 1

        print("numver_of_element: ", len(Elements))
        for element in Elements:
            # Collect valid integration points
            integration_points = [
                [point.X(), point.Y(), point.Z(), point.Weight()]
                for point in element.GetIntegrationPoints()
                if(point.Weight() > 0)
            ]

            # Create elements
            if integration_points:
                el_count += 1
                # Create quadrature_point_geometries
                quadrature_point_geometries = KM.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2, integration_points)
                KratosNurbsVolumeModelPart.CreateNewElement(
                    element_name,
                    el_count,
                    quadrature_point_geometries[0],
                    volume_properties
                )

        return el_count

    @staticmethod
    def add_conditions_to_model_part(
            KratosNurbsVolumeModelPart: KM.ModelPart,
            Conditions,
            BoundsXYZ: Tuple[Point3D, Point3D],
            BoundsUVW: Tuple[Point3D, Point3D],
            nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        """
        Adds the QuESo conditions to the Kratos NurbsVolume ModelPart.

        Args:
            KratosNurbsVolumeModelPart (KM.ModelPart): The Kratos model part to which the conditions will be added.
            Conditions (QuESo.ConditionVector): List of QuESo boundary conditions to add to the model part.
            BoundsXYZ (Tuple[Point3D, Point3D]): Lower and upper bounds in XYZ coordinates for the conditions.
            BoundsUVW (Tuple[Point3D, Point3D]): Lower and upper bounds in UVW coordinates for the conditions.
        """

        boundary_conditions = []
        for bc in Conditions:
            if isinstance(bc, dict):
                if "segments" not in bc:
                    raise RuntimeError("Boundary condition definition must contain 'segments'.")
                for segment in bc["segments"]:
                    weak_bc = ModelPartUtilities._create_bc_from_settings(bc, segment, nurbs_volume_name)
                    weak_bc.bounds_xyz = BoundsXYZ
                    weak_bc.bounds_uvw = BoundsUVW
                    boundary_conditions.append(weak_bc)
            elif bc.is_weak_condition():
                condition_settings = bc.GetSettings()
                type_name = condition_settings.GetString("condition_type")

                for segment in bc:
                    weak_bc = ModelPartUtilities._create_bc_from_settings(
                        {
                            "condition_type": type_name,
                            "value": condition_settings.GetDoubleVector("value") if condition_settings.Has("value") else [0.0, 0.0, 0.0],
                            "penalty_factor": condition_settings.GetDouble("penalty_factor") if condition_settings.Has("penalty_factor") else 0.0,
                            "modulus": condition_settings.GetDouble("modulus") if condition_settings.Has("modulus") else 0.0,
                            "direction": condition_settings.GetDoubleVector("direction") if condition_settings.Has("direction") else [0.0, 0.0, 0.0],
                        },
                        segment,
                        nurbs_volume_name,
                    )
                    weak_bc.bounds_xyz = BoundsXYZ
                    weak_bc.bounds_uvw = BoundsUVW
                    boundary_conditions.append(weak_bc)
            else:
                boundary_conditions.append(bc)

        for bc in boundary_conditions:
            bc.apply(KratosNurbsVolumeModelPart)

    @staticmethod
    def remove_all_elements(KratosModelPart: KM.ModelPart) -> None:
        """
        Removes all elements from the Kratos ModelPart.

        Args:
            KratosModelPart (KM.ModelPart): The Kratos ModelPart from which all elements will be removed.
        """
        for element in KratosModelPart.Elements:
            element.Set(KM.TO_ERASE, True)
        KratosModelPart.RemoveElements(KM.TO_ERASE)

    @staticmethod
    def remove_all_conditions(KratosModelPart: KM.ModelPart) -> None:
        """
        Removes all conditions from the Kratos ModelPart.

        Args:
            KratosModelPart (KM.ModelPart): The Kratos ModelPart from which all conditions will be removed.
        """
        for condition in KratosModelPart.Conditions:
            condition.Set(KM.TO_ERASE, True)
        KratosModelPart.RemoveConditions(KM.TO_ERASE)
