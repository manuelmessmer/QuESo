import KratosMultiphysics
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
        ''' Reads Kratos ModelPart from the TIBRA triangle mesh. '''
        vertices = TriangleMesh.GetVertices()
        for v_id, vertex in enumerate(vertices):
            KratosModelPart.CreateNewNode(v_id+1, vertex[0], vertex[1], vertex[2])

        triangles = TriangleMesh.GetTriangles()
        for t_id, triangle in enumerate(triangles):
            KratosModelPart.CreateNewElement("ShellThinElement3D3N", t_id+1, [triangle[0]+1, triangle[1]+1, triangle[2]+1], KratosModelPart.GetProperties()[1])

    @staticmethod
    def CreateTIBRAInput(KratosEmbeddedModelPart, TibraParameters):
        ''' Writes the KratosEmbeddedModelPart (including submodelpart for conditions) to STL files, which can be read by TIRBA. '''

        # Write main model part
        input_filename = TibraParameters.GetInputFilename()
        m_lower_case = re.search('/(.*).stl', input_filename)
        m_upper_case = re.search('/(.*).STL', input_filename)
        if m_lower_case:
            model_part_name = m_lower_case.group(1)
        elif m_upper_case:
            model_part_name = m_lower_case.group(1)
        else:
            raise Exception("CreateTIBRAInput::Filename is not valid.")

        ModelPartUtilities._WriteModelPartToSTL(KratosEmbeddedModelPart.GetSubModelPart(model_part_name), input_filename)

        # Write condition model part
        condition_filenames = {}
        for cond_id in range(TibraParameters.NumberOfConditions()):
            tmp_filename = TibraParameters.GetFilenameOfCondition(cond_id)
            m_lower_case = re.search('/(.*).stl', tmp_filename)
            m_upper_case = re.search('/(.*).STL', tmp_filename)
            if m_lower_case:
                condition_filenames[m_lower_case.group(1)] = cond_id
            elif m_upper_case:
                condition_filenames[m_lower_case.group(1)] = cond_id
            else:
                raise Exception("CreateTIBRAInput::Filename is not valid.")

        for sub_model_part in KratosEmbeddedModelPart.SubModelParts:
            sub_model_part_name = sub_model_part.Name
            if( sub_model_part_name in condition_filenames.keys()):
                cond_id = condition_filenames[sub_model_part_name]
                ModelPartUtilities._WriteModelPartToSTL(sub_model_part, TibraParameters.GetFilenameOfCondition(cond_id))

    @staticmethod
    def AddElementsToModelPart(KratosNurbsVolumeModelPart, Elements):
        ''' Adds the TIBRA elements to the KratosNurbsVolumeModelPart. '''
        nurbs_volume = KratosNurbsVolumeModelPart.GetGeometry("NurbsVolume")
        volume_properties = KratosNurbsVolumeModelPart.GetProperties()[1]

        el_count = 0
        for element in Elements:
            integration_points = []
            if element.IsTrimmed():
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
                KratosNurbsVolumeModelPart.CreateNewElement('SmallDisplacementElement3D8N', el_count, quadrature_point_geometries[0], volume_properties)

    @staticmethod
    def AddConditionsToModelPart(KratosNurbsVolumeModelPart, Conditions, LowerBound, UpperBound):
        ''' Adds the TIBRA elements to the KratosNurbsVolumeModelPart. '''
        boundary_conditions = []
        for bc in Conditions:
            if bc.IsWeakCondition():
                if( bc.Type() == "dirichlet" ):
                    dirichlet_triangles = bc.GetTriangleMesh()
                    boundary_conditions.append(
                        PenaltySupport(dirichlet_triangles, LowerBound, UpperBound, bc.GetPrescribed(), bc.GetPenaltyFactor()) )
                elif( bc.Type() == "neumann" ):
                    neumann_triangles = bc.GetTriangleMesh()
                    boundary_conditions.append(
                        SurfaceLoad(neumann_triangles, LowerBound, UpperBound, bc.GetPrescribed(), False) )
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
