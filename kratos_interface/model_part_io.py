
import KratosMultiphysics as KM

def WriteKratosModelPart(TibraTriangleMesh, KratosModelPart):
    vertices = TibraTriangleMesh.GetVertices()
    for v_id, vertex in enumerate(vertices):
        KratosModelPart.CreateNewNode(v_id+1, vertex[0], vertex[1], vertex[2])

    triangles = TibraTriangleMesh.GetTriangles()
    for t_id, triangle in enumerate(triangles):
        KratosModelPart.CreateNewElement("ShellThinElement3D3N", t_id+1, [triangle[0]+1, triangle[1]+1, triangle[2]+1], KratosModelPart.GetProperties()[1])

def ReadKratosModelPart(KratosModelPart, TibraTriangleMesh):
    for node in KratosModelPart.Nodes:
        node.Id()