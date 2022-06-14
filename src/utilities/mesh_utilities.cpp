#include <vtkSTLReader.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataNormals.h>


#include <string>
#include <iostream>
#include "includes/isotropic_remeshing/isotropicremesher.h"
#include "includes/isotropic_remeshing/isotropichalfedgemesh.h"
#include "includes/isotropic_remeshing/vector3.h"

#include "includes/vtkbool/vtkPolyDataBooleanFilter.h"
#include "includes/vtkbool/vtkPolyDataContactFilter.h"

#include "utilities/mesh_utilities.h"
#include "io/io_utilities.h"


std::chrono::duration<double> Elapsed_Time{};
int Point::_tag = 0; // Required for vtkbool

bool MeshUtilities::DoIntersect(vtkSmartPointer<vtkPolyData> pInput1, vtkSmartPointer<vtkPolyData> pInput2){
    vtkSmartPointer<vtkPolyDataContactFilter> cl = vtkSmartPointer<vtkPolyDataContactFilter>::New();

    cl->SetInputData(0, pInput1);
    cl->SetInputData(1, pInput2);
    cl->Update();

    if (cl->GetOutput()->GetNumberOfCells() == 0) {
        cl->RemoveAllInputs();
        return 0;
    }
    cl->RemoveAllInputs();
    return 1;
}

vtkSmartPointer<vtkPolyData> MeshUtilities::IsotropicRemeshing(vtkSmartPointer<vtkPolyData> pInput, const Parameters& rParam, int id){

        vtkSmartPointer<vtkTriangleFilter> triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();

        triangle_filter->SetInputDataObject(pInput);
        triangle_filter->Update();
        auto intersection_mesh_tris = triangle_filter->GetOutput();

        vtkSmartPointer<vtkPoints> points = intersection_mesh_tris->GetPoints();
        vtkSmartPointer<vtkDataArray> dataArray = points->GetData();
        const int number_of_points = intersection_mesh_tris->GetNumberOfPoints();
        const int number_of_triangles = intersection_mesh_tris->GetNumberOfPolys();

        std::vector<Vector3> inputVertices{};
        inputVertices.reserve(number_of_points);
        std::vector<std::vector<size_t>> inputTriangles{};
        inputTriangles.reserve(number_of_triangles);
        for( int i = 0; i < number_of_points; i++){
            inputVertices.push_back( Vector3{dataArray->GetComponent(i,0), dataArray->GetComponent(i,1), dataArray->GetComponent(i,2)} );
        }

        auto polys = intersection_mesh_tris->GetPolys();
        polys->InitTraversal();
        for( int i = 0; i < number_of_triangles; i++){
            vtkIdType npts, *pts;
            auto h = polys->GetNextCell(npts, pts);
            if( h == 0 ){
                break;
            }
            if( npts == 3 ){
                std::vector<size_t> tmp_vec = { (size_t)pts[0], (size_t)pts[1], (size_t)pts[2]};
                inputTriangles.push_back( tmp_vec );
            } else {
                std::cout << "MeshUtilities :: npts != 3" << std::endl;
            }
        }

        // Initialize Remesher
        double value = 1.0;

        auto start_time = std::chrono::high_resolution_clock::now();
        IsotropicRemesher isotropicRemesher(&inputVertices, &inputTriangles);
        isotropicRemesher.setSharpEdgeIncludedAngle(91);
        isotropicRemesher.setTargetEdgeLength(isotropicRemesher.initialAverageEdgeLength() * value);
        //isotropicRemesher.initialAverageEdgeLength
        //isotropicRemesher.setTargetTriangleCount(rParam.MinimumNumberOfTriangles());
        isotropicRemesher.remesh(5);

        int poly_count = 0;
        IsotropicHalfedgeMesh *halfedgeMesh = isotropicRemesher.remeshedHalfedgeMesh();
        for (IsotropicHalfedgeMesh::Face *face = halfedgeMesh->moveToNextFace(nullptr);
                    nullptr != face;
                    face = halfedgeMesh->moveToNextFace(face)) {
            poly_count++;
        }

        int iteration = 0;
        while( poly_count < rParam.MinimumNumberOfTriangles() && iteration < 10) {
            value *= 0.95 * std::sqrt( (double)poly_count / (double)rParam.MinimumNumberOfTriangles() );

            isotropicRemesher.setTargetEdgeLength(isotropicRemesher.initialAverageEdgeLength() * value );
            isotropicRemesher.remesh(2);
            halfedgeMesh = isotropicRemesher.remeshedHalfedgeMesh();
            poly_count = 0;
            for (IsotropicHalfedgeMesh::Face *face = halfedgeMesh->moveToNextFace(nullptr);
                    nullptr != face;
                    face = halfedgeMesh->moveToNextFace(face)) {
                poly_count++;
            }
            //std::cout << "Ploy count: " << poly_count << std::endl;
            iteration++;
        }


        auto end_time = std::chrono::high_resolution_clock::now();
        Elapsed_Time += end_time - start_time;

        size_t outputIndex = 0;

        vtkSmartPointer<vtkPoints> new_points = vtkSmartPointer<vtkPoints>::New();
        for (IsotropicHalfedgeMesh::Vertex *vertex = halfedgeMesh->moveToNextVertex(nullptr);
                nullptr != vertex;
                vertex = halfedgeMesh->moveToNextVertex(vertex)) {

            vertex->outputIndex = outputIndex++;
            new_points->InsertNextPoint( vertex->position[0], vertex->position[1], vertex->position[2] );
        }

        vtkSmartPointer<vtkCellArray> new_triangles = vtkSmartPointer<vtkCellArray>::New();
        for (IsotropicHalfedgeMesh::Face *face = halfedgeMesh->moveToNextFace(nullptr);
                nullptr != face;
                face = halfedgeMesh->moveToNextFace(face)) {

            vtkSmartPointer<vtkTriangle> tmp_triangle = vtkSmartPointer<vtkTriangle>::New();
            tmp_triangle->GetPointIds()->SetId(0, (size_t)face->halfedge->previousHalfedge->startVertex->outputIndex );
            tmp_triangle->GetPointIds()->SetId(1, (size_t)face->halfedge->startVertex->outputIndex);
            tmp_triangle->GetPointIds()->SetId(2, (size_t)face->halfedge->nextHalfedge->startVertex->outputIndex);

            new_triangles->InsertNextCell( tmp_triangle );
        }

        vtkSmartPointer<vtkPolyData> remeshed = vtkSmartPointer<vtkPolyData>::New();
        remeshed->SetPoints( new_points );
        remeshed->SetPolys( new_triangles );


        // std::string str = "remeshed";
        // str.append( std::to_string(id) );
        //str.append(".vtk");
        //IO_Utilities::WriteVTK(remeshed, "output.vtk");
        return remeshed;
}

vtkSmartPointer<vtkPolyData> MeshUtilities::GetIntersection(vtkSmartPointer<vtkPolyData> pInput1, vtkSmartPointer<vtkPolyData> pInput2, const Parameters& rParam){

        // Get vtkBool filter and compute intersection mesh
        //std::cout << "Here... " << std::endl;
        vtkSmartPointer<vtkPolyDataBooleanFilter> bool_filter = vtkPolyDataBooleanFilter::New();
        vtkSmartPointer<vtkPolyData> pInput1_copy = vtkSmartPointer<vtkPolyData>::New();
        pInput1_copy->DeepCopy(pInput1);

        bool_filter->SetInputData(0, pInput1_copy);
        bool_filter->SetInputData(1, pInput2);
        bool_filter->SetOperModeToIntersection();
        bool_filter->Update();

        vtkSmartPointer<vtkPolyData> intersection_mesh = bool_filter->GetOutput();
        //std::cout << "After intersection! " << std::endl;

        //std::cout << "Here...22 " << std::endl;
        //IO_Utilities::WriteVTK(intersection_mesh, "inter.vtk");
        return intersection_mesh;
    }