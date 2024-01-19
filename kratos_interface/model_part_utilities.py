import KratosMultiphysics
import KratosMultiphysics.OptimizationApplication as OptimizationApplication
from kratos_interface.weak_bcs import PenaltySupport
from kratos_interface.weak_bcs import SurfaceLoad
import re

class ModelPartUtilities:
    @staticmethod
    def _WriteModelPartToSTL(KratosModelPart, Filename):
        ''' Writes KratosModelPart to STL. '''
        io_settings = KratosMultiphysics.Parameters(
            """{
                "open_mode" : "write",
                "new_entity_type" : "geometry"
            }""")

        stl_io = KratosMultiphysics.StlIO(Filename, io_settings)
        stl_io.WriteModelPart(KratosModelPart)

    @staticmethod
    def ReadModelPartFromTriangleMesh(KratosModelPart, TriangleMesh):
        ''' Reads Kratos ModelPart from the QuESo triangle mesh. '''
        vertices = TriangleMesh.GetVertices()
        for v_id, vertex in enumerate(vertices):
            KratosModelPart.CreateNewNode(v_id+1, vertex[0], vertex[1], vertex[2])

        triangles = TriangleMesh.GetTriangles()
        for t_id, triangle in enumerate(triangles):
            KratosModelPart.CreateNewElement("ShellThinElement3D3N", t_id+1, [triangle[0]+1, triangle[1]+1, triangle[2]+1], KratosModelPart.GetProperties()[1])

    @staticmethod
    def CreateQuESoInput(KratosEmbeddedModelPart, QuESoParameters):
        ''' Writes the KratosEmbeddedModelPart (including submodelpart for conditions) to STL files, which can be read by TIRBA. '''

        # Write main model part
        input_filename = QuESoParameters.GetInputFilename()
        m_lower_case = re.search('/(.*).stl', input_filename)
        m_upper_case = re.search('/(.*).STL', input_filename)
        if m_lower_case:
            model_part_name = m_lower_case.group(1)
        elif m_upper_case:
            model_part_name = m_lower_case.group(1)
        else:
            raise Exception("CreateQuESoInput::Filename is not valid.")

        ModelPartUtilities._WriteModelPartToSTL(KratosEmbeddedModelPart.GetSubModelPart(model_part_name), input_filename)

        # Write condition model part
        condition_filenames = {}
        for cond_id in range(QuESoParameters.NumberOfConditions()):
            tmp_filename = QuESoParameters.GetFilenameOfCondition(cond_id)
            m_lower_case = re.search('/(.*).stl', tmp_filename)
            m_upper_case = re.search('/(.*).STL', tmp_filename)
            if m_lower_case:
                condition_filenames[m_lower_case.group(1)] = cond_id
            elif m_upper_case:
                condition_filenames[m_lower_case.group(1)] = cond_id
            else:
                raise Exception("CreateQuESoInput::Filename is not valid.")

        for sub_model_part in KratosEmbeddedModelPart.SubModelParts:
            sub_model_part_name = sub_model_part.Name
            if( sub_model_part_name in condition_filenames.keys()):
                cond_id = condition_filenames[sub_model_part_name]
                ModelPartUtilities._WriteModelPartToSTL(sub_model_part, QuESoParameters.GetFilenameOfCondition(cond_id))

    @staticmethod
    def AddElementsToModelPart(KratosNurbsVolumeModelPart, Elements):
        ''' Adds the QuESo elements to the KratosNurbsVolumeModelPart. '''
        nurbs_volume = KratosNurbsVolumeModelPart.GetGeometry("NurbsVolume")
        volume_properties = KratosNurbsVolumeModelPart.GetProperties()[1]
        import QuESo_PythonApplication as QuESoApp
        from QuESo_PythonApplication.PyQuESo import PyQuESo

        pyqueso = PyQuESo("QUESOParameters_tmp.json")
        nodes = QuESoApp.PointVector()
        for node in KratosNurbsVolumeModelPart.Nodes:
            nodes.append( QuESoApp.Point(node.X0, node.Y0, node.Z0) )

        pyqueso.Run()
        are_inside = pyqueso.IsInside(nodes)
        print("Number of nodes: ---------------------------------------------")
        print("Number of nodes: ", len(nurbs_volume))
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.TO_ERASE, False, KratosNurbsVolumeModelPart.Nodes)
        for node, is_inside in zip(KratosNurbsVolumeModelPart.Nodes, are_inside):
            if( is_inside ):
                node.Set(KratosMultiphysics.TO_ERASE, True)

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.TO_ERASE, False, KratosNurbsVolumeModelPart.Elements)
        el_count = 0
        for element in Elements:
            integration_points = []
            delta_low = element.LowerBoundUVW()
            delta_up = element.UpperBoundUVW()
            upper_bound = element.UpperBoundXYZ()
            lower_bound = element.LowerBoundXYZ()
            delta_x = upper_bound[0] - lower_bound[0]
            delta_y = upper_bound[1] - lower_bound[1]
            delta_z = upper_bound[2] - lower_bound[2]
            #extend = [0.0, 0.0, 0.0]
            trimmed = False
            if element.IsTrimmed():
                trimmed = True
                for point in element.GetIntegrationPoints():
                    weight = point.GetWeight()
                    if( weight > 0):
                        integration_points.append([point.GetX(), point.GetY(), point.GetZ(), point.GetWeight()])
            else:
                for point in element.GetIntegrationPoints():
                    integration_points.append([point.GetX(), point.GetY(), point.GetZ(), point.GetWeight()])

            if( len(integration_points) > 0 ):
                el_count += 1
                # Create quadrature_point_geometries
                quadrature_point_geometries = KratosMultiphysics.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2, integration_points)
                kratos_element = KratosNurbsVolumeModelPart.CreateNewElement('SmallDisplacementElement3D8N', el_count, quadrature_point_geometries[0], volume_properties)
                kratos_element.Set(KratosMultiphysics.TO_ERASE, trimmed)
                import numpy as np
                num_points = 12
                if trimmed:
                    triangle_mesh = element.GetClippedTriangleMesh()
                    num_of_triangles = triangle_mesh.NumOfTriangles()
                    integration_points = []

                    matrix = KratosMultiphysics.Matrix(num_of_triangles*num_points, 7)
                    for triangle_id in range(num_of_triangles):

                        points = triangle_mesh.GetIntegrationPointsGlobal(triangle_id, 3)
                        normal = triangle_mesh.Normal(triangle_id)
                        # if( triangle_id == 0 ):
                        #     print(points)
                        #     #print(points.GetWeight())
                        #     ghdgf55

                        for point_id, point in enumerate(points):
                            matrix[num_points*triangle_id+point_id, 0] = point.GetX()
                            matrix[num_points*triangle_id+point_id, 1] = point.GetY()
                            matrix[num_points*triangle_id+point_id, 2] = point.GetZ()
                            matrix[num_points*triangle_id+point_id, 3] = point.GetWeight()
                            matrix[num_points*triangle_id+point_id, 4] = normal[0]
                            matrix[num_points*triangle_id+point_id, 5] = normal[1]
                            matrix[num_points*triangle_id+point_id, 6] = normal[2]

                    kratos_element.SetValue(OptimizationApplication.EMBEDDED_MESH,  matrix)
                point1 = KratosMultiphysics.Vector([0.0, 0.0, 0.0])
                point1[0] = delta_x
                point1[1] = delta_y
                point1[2] = delta_z
                kratos_element.SetValue(KratosMultiphysics.CONTACT_FORCE, point1)

                point2 = KratosMultiphysics.Vector([0.0, 0.0, 0.0])
                point2[0] = delta_low[0]
                point2[1] = delta_low[1]
                point2[2] = delta_low[2]
                kratos_element.SetValue(KratosMultiphysics.INTERNAL_FORCE, point2)

                point3 = KratosMultiphysics.Vector([0.0, 0.0, 0.0])
                point3[0] = delta_up[0]
                point3[1] = delta_up[1]
                point3[2] = delta_up[2]
                kratos_element.SetValue(KratosMultiphysics.EXTERNAL_FORCE, point3)

    @staticmethod
    def AddConditionsToModelPart(KratosNurbsVolumeModelPart, Conditions, BoundsXYZ, BoundsUVW):
        ''' Adds the QuESo elements to the KratosNurbsVolumeModelPart. '''
        boundary_conditions = []
        for bc in Conditions:
            if bc.IsWeakCondition():
                if( bc.Type() == "dirichlet" ):
                    dirichlet_triangles = bc.GetTriangleMesh()
                    boundary_conditions.append(
                        PenaltySupport(dirichlet_triangles, BoundsXYZ, BoundsUVW, bc.GetPrescribed(), bc.GetPenaltyFactor()) )
                elif( bc.Type() == "neumann" ):
                    neumann_triangles = bc.GetTriangleMesh()
                    boundary_conditions.append(
                        SurfaceLoad(neumann_triangles, BoundsXYZ, BoundsUVW, bc.GetPrescribed(), False) )
            else:
                boundary_conditions.append(bc)

        for bc in boundary_conditions:
            bc.apply(KratosNurbsVolumeModelPart)

    @staticmethod
    def RemoveAllElements(KratosModelPart):
        ''' Removes all elements from the KratosModelPart. '''
        for element in KratosModelPart.Elements:
            element.Set(KratosMultiphysics.TO_ERASE, True)
        KratosModelPart.RemoveElements(KratosMultiphysics.TO_ERASE)

    @staticmethod
    def RemoveAllConditions(KratosModelPart):
        ''' Removes all conditions from the KratosModelPart. '''
        for condition in KratosModelPart.Conditions:
            condition.Set(KratosMultiphysics.TO_ERASE, True)
        KratosModelPart.RemoveConditions(KratosMultiphysics.TO_ERASE)
