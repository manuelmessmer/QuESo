from typing import List, Literal
# Import Kratos
import KratosMultiphysics as KM
# Import QuESo modules
import QuESo_PythonApplication as QuESo
from kratos_interface.weak_bcs import PenaltySupport
from kratos_interface.weak_bcs import LagrangeSupport
from kratos_interface.weak_bcs import SurfaceLoad
from kratos_interface.weak_bcs import PressureLoad

class ModelPartUtilities:
    """
    Utility class for managing Kratos ModelPart operations related to QuESo triangle meshes and elements.

    This class provides various static methods to interact with Kratos model parts, including reading from and writing
    to STL files, transferring triangle meshes between QuESo and Kratos, and adding elements and conditions to Kratos
    model parts.
    """

    @staticmethod
    def _WriteModelPartToSTL(
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

        stl_io = KM.StlIO(Filename, io_settings)
        stl_io.WriteModelPart(KratosModelPart)

    @staticmethod
    def ReadModelPartFromTriangleMesh(
            KratosModelPart: KM.ModelPart,
            TriangleMesh: QuESo.TriangleMesh
        ) -> None:
        """
        Reads a Kratos ModelPart from a QuESo triangle mesh.

        Args:
            KratosModelPart (KM.ModelPart): The Kratos ModelPart to populate with data from the triangle mesh.
            TriangleMesh (TriangleMesh): The QuESo triangle mesh from which the model part is read.
        """
        vertices = TriangleMesh.GetVertices()
        for v_id, vertex in enumerate(vertices):
            KratosModelPart.CreateNewNode(v_id+1, vertex[0], vertex[1], vertex[2])

        triangles = TriangleMesh.GetTriangles()
        for t_id, triangle in enumerate(triangles):
            KratosModelPart.CreateNewElement("ShellThinElement3D3N", t_id+1, [triangle[0]+1, triangle[1]+1, triangle[2]+1], KratosModelPart.GetProperties()[1])

    @staticmethod
    def ReadTriangleMeshFromModelPart(
            TriangleMesh: QuESo.TriangleMesh,
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
        for node in KratosModelPart.Nodes:
            TriangleMesh.AddVertex([node.X, node.Y, node.Z])

        if( type == "Elements"):
            entity_list = KratosModelPart.Elements
        elif(type == "Conditions"):
            entity_list = KratosModelPart.Conditions
        else:
            message = "ModelPartUtilities :: ReadTriangleMeshFromModelPart :: Given type: '" + str(type)
            message += "' not valid. Available options are: 'Elements' and 'Conditions'."
            raise Exception(message)

        for entity in entity_list:
            geometry = entity.GetGeometry()

            # Ensure the geometry is a triangle
            if( len(geometry) != 3 ):
                raise Exception("ModelPartUtilities :: ReadTriangleMeshFromModelPart :: Queso only allows triangles.")

            # Get the node ids for the triangle
            node_ids = [id_map[node.Id] for node in geometry]
            TriangleMesh.AddTriangle(node_ids)

            # Compute the normal for the triangle
            normal = QuESo.Triangle.NormalStatic(
                [node.X, node.Y, node.Z] for node in geometry
            )
            TriangleMesh.AddNormal(normal)

    @staticmethod
    def AddElementsToModelPart(
            KratosNurbsVolumeModelPart: KM.ModelPart,
            Elements: QuESo.ElementVector
        ) -> None:
        """
        Adds the QuESo elements to the Kratos NurbsVolume ModelPart.

        Args:
            KratosNurbsVolumeModelPart (KM.ModelPart): The Kratos model part to which the elements will be added.
            Elements (QuESo.ElementVector): List of QuESo elements to add to the model part.
        """
        nurbs_volume = KratosNurbsVolumeModelPart.GetGeometry("NurbsVolume")
        volume_properties = KratosNurbsVolumeModelPart.GetProperties()[1]

        el_count = 0
        for element in Elements:
            integration_points = []
            if element.IsTrimmed():
                for point in element.GetIntegrationPoints():
                    weight = point.Weight()
                    if( weight > 0):
                        integration_points.append([point.X(), point.Y(), point.Z(), point.Weight()])
            else:
                for point in element.GetIntegrationPoints():
                    integration_points.append([point.X(), point.Y(), point.Z(), point.Weight()])

            if( len(integration_points) > 0 ):
                el_count += 1
                # Create quadrature_point_geometries
                quadrature_point_geometries = KM.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2, integration_points)
                KratosNurbsVolumeModelPart.CreateNewElement('SmallDisplacementElement3D8N', el_count, quadrature_point_geometries[0], volume_properties)

    @staticmethod
    def AddConditionsToModelPart(
            KratosNurbsVolumeModelPart: KM.ModelPart,
            Conditions: QuESo.ConditionVector,
            BoundsXYZ: List[List[float]],
            BoundsUVW: List[List[float]]
        ) -> None:
        """
        Adds the QuESo conditions to the Kratos NurbsVolume ModelPart.

        Args:
            KratosNurbsVolumeModelPart (KM.ModelPart): The Kratos model part to which the conditions will be added.
            Conditions (QuESo.ConditionVector): List of QuESo boundary conditions to add to the model part.
            BoundsXYZ (List[List[float]]): Lower and upper bounds in XYZ coordinates for the conditions.
            BoundsUVW (List[List[float]]): Lower and upper bounds in UVW coordinates for the conditions.
        """

        # Define a mapping of condition types to their handlers
        condition_handlers = {
            "PenaltySupportCondition": lambda settings, segment: PenaltySupport(
                segment.GetTriangleMesh(),
                BoundsXYZ,
                BoundsUVW,
                settings.GetDoubleVector("value"),
                settings.GetDouble("penalty_factor")
            ),
            "LagrangeSupportCondition": lambda settings, segment: LagrangeSupport(
                segment.GetTriangleMesh(),
                BoundsXYZ,
                BoundsUVW,
                settings.GetDoubleVector("value")
            ),
            "SurfaceLoadCondition": lambda settings, segment: SurfaceLoad(
                segment.GetTriangleMesh(),
                BoundsXYZ,
                BoundsUVW,
                settings.GetDouble("modulus"),
                settings.GetDoubleVector("direction")
            ),
            "PressureLoadCondition": lambda settings, segment: PressureLoad(
                segment.GetTriangleMesh(),
                BoundsXYZ,
                BoundsUVW,
                settings.GetDouble("modulus")
            )
        }

        boundary_conditions = []
        for bc in Conditions:
            if bc.IsWeakCondition():
                condition_settings = bc.GetSettings()
                type_name = condition_settings.GetString("condition_type")

                if type_name not in condition_handlers:
                    raise Exception(f"Given condition type '{type_name}' is not available.")

                handler = condition_handlers[type_name]
                for segment in bc:
                    boundary_conditions.append(handler(condition_settings, segment))
            else:
                boundary_conditions.append(bc)

        for bc in boundary_conditions:
            bc.apply(KratosNurbsVolumeModelPart)

    @staticmethod
    def RemoveAllElements(KratosModelPart: KM.ModelPart) -> None:
        """
        Removes all elements from the Kratos ModelPart.

        Args:
            KratosModelPart (KM.ModelPart): The Kratos ModelPart from which all elements will be removed.
        """
        for element in KratosModelPart.Elements:
            element.Set(KM.TO_ERASE, True)
        KratosModelPart.RemoveElements(KM.TO_ERASE)

    @staticmethod
    def RemoveAllConditions(KratosModelPart: KM.ModelPart) -> None:
        """
        Removes all conditions from the Kratos ModelPart.

        Args:
            KratosModelPart (KM.ModelPart): The Kratos ModelPart from which all conditions will be removed.
        """
        for condition in KratosModelPart.Conditions:
            condition.Set(KM.TO_ERASE, True)
        KratosModelPart.RemoveConditions(KM.TO_ERASE)
