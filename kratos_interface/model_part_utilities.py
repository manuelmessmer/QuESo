
# Import QuESo
import QuESo_PythonApplication as QuESo_App
# Import QuESo
import KratosMultiphysics
from kratos_interface.weak_bcs import PenaltySupport
from kratos_interface.weak_bcs import LagrangeSupport
from kratos_interface.weak_bcs import SurfaceLoad
from kratos_interface.weak_bcs import PressureLoad
from QuESo_PythonApplication import TriangleMesh as TriangleUtilities

class ModelPartUtilities:
    @staticmethod
    def _WriteModelPartToSTL(KratosModelPart, Filename):
        ''' Deprectaed method.
        Writes KratosModelPart to STL.
        '''
        io_settings = KratosMultiphysics.Parameters(
            """{
                "open_mode" : "write",
                "new_entity_type" : "geometry"
            }""")

        stl_io = KratosMultiphysics.StlIO(Filename, io_settings)
        stl_io.WriteModelPart(KratosModelPart)

    @staticmethod
    def ReadModelPartFromTriangleMesh(KratosModelPart, TriangleMesh):
        ''' Reads Kratos ModelPart from the QuESo triangle mesh.
        '''
        vertices = TriangleMesh.GetVertices()
        for v_id, vertex in enumerate(vertices):
            KratosModelPart.CreateNewNode(v_id+1, vertex[0], vertex[1], vertex[2])

        triangles = TriangleMesh.GetTriangles()
        for t_id, triangle in enumerate(triangles):
            KratosModelPart.CreateNewElement("ShellThinElement3D3N", t_id+1, [triangle[0]+1, triangle[1]+1, triangle[2]+1], KratosModelPart.GetProperties()[1])

    @staticmethod
    def ReadTriangleMeshFromModelPart(TriangleMesh, KratosModelPart, type="Elements"):
        ''' Reads QuESo triangle mesh from Kratos ModelPart
        '''
        id_map = {}
        for queso_id, node in enumerate(KratosModelPart.Nodes):
            kratos_id = node.Id
            id_map[kratos_id] = queso_id
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
            if( len(geometry) != 3 ):
                raise Exception("ModelPartUtilities :: ReadTriangleMeshFromModelPart :: Queso only allows triangles.")
            node1 = geometry[0]
            node2 = geometry[1]
            node3 = geometry[2]
            triangle = [id_map[node1.Id], id_map[node2.Id], id_map[node3.Id] ]
            TriangleMesh.AddTriangle(triangle)
            normal = TriangleUtilities.NormalStatic( [node1.X, node1.Y, node1.Z], [node2.X, node2.Y, node2.Z], [node3.X, node3.Y, node3.Z] )
            TriangleMesh.AddNormal(normal)

    @staticmethod
    def AddElementsToModelPart(KratosNurbsVolumeModelPart, Elements):
        ''' Adds the QuESo elements to the KratosNurbsVolumeModelPart.
        '''
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
                quadrature_point_geometries = KratosMultiphysics.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2, integration_points)
                KratosNurbsVolumeModelPart.CreateNewElement('SmallDisplacementElement3D8N', el_count, quadrature_point_geometries[0], volume_properties)

    @staticmethod
    def AddConditionsToModelPart(KratosNurbsVolumeModelPart, Conditions, BoundsXYZ, BoundsUVW):
        ''' Adds the QuESo elements to the KratosNurbsVolumeModelPart. '''
        boundary_conditions = []
        for bc in Conditions:
            if bc.IsWeakCondition():
                condition_settings = bc.GetSettings()
                type_name = condition_settings.GetString(QuESo_App.ConditionSettings.condition_type)
                if( type_name == "PenaltySupportCondition" ):
                    prescribed_displacement = condition_settings.GetDoubleVector(QuESo_App.ConditionSettings.value)
                    penalty_factor = condition_settings.GetDouble(QuESo_App.ConditionSettings.penalty_factor)
                    for condition_segment in bc:
                        dirichlet_triangles = condition_segment.GetTriangleMesh()
                        boundary_conditions.append(PenaltySupport(dirichlet_triangles, BoundsXYZ, BoundsUVW, prescribed_displacement, penalty_factor) )
                elif( type_name == "LagrangeSupportCondition" ):
                    prescribed_displacement = condition_settings.GetDoubleVector(QuESo_App.ConditionSettings.value)
                    for condition_segment in bc:
                        dirichlet_triangles = condition_segment.GetTriangleMesh()
                        boundary_conditions.append(LagrangeSupport(dirichlet_triangles, BoundsXYZ, BoundsUVW, prescribed_displacement) )
                elif( type_name == "SurfaceLoadCondition" ):
                    modulus = condition_settings.GetDouble(QuESo_App.ConditionSettings.modulus)
                    direction = condition_settings.GetDoubleVector(QuESo_App.ConditionSettings.direction)
                    for condition_segment in bc:
                        neumann_triangles = condition_segment.GetTriangleMesh()
                        boundary_conditions.append(SurfaceLoad(neumann_triangles, BoundsXYZ, BoundsUVW, modulus, direction) )
                elif( type_name == "PressureLoadCondition" ):
                    modulus = condition_settings.GetDouble(QuESo_App.ConditionSettings.modulus)
                    for condition_segment in bc:
                        neumann_triangles = condition_segment.GetTriangleMesh()
                        boundary_conditions.append(PressureLoad(neumann_triangles, BoundsXYZ, BoundsUVW, modulus) )
                else:
                    message = "Given condition type '" + type_name + "' is not available.\n"
                    raise Exception(message)
            else:
                boundary_conditions.append(bc)

        for bc in boundary_conditions:
            bc.apply(KratosNurbsVolumeModelPart)

    @staticmethod
    def RemoveAllElements(KratosModelPart):
        ''' Removes all elements from the KratosModelPart.
        '''
        for element in KratosModelPart.Elements:
            element.Set(KratosMultiphysics.TO_ERASE, True)
        KratosModelPart.RemoveElements(KratosMultiphysics.TO_ERASE)

    @staticmethod
    def RemoveAllConditions(KratosModelPart):
        ''' Removes all conditions from the KratosModelPart.
        '''
        for condition in KratosModelPart.Conditions:
            condition.Set(KratosMultiphysics.TO_ERASE, True)
        KratosModelPart.RemoveConditions(KratosMultiphysics.TO_ERASE)
