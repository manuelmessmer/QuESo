// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// VTK includes
#include <vtkCubeSource.h>


//// Project includes
#include "TrIGA_main.hpp"
#include "utilities/integration_point_utilities.h"
#include "utilities/moment_fitting_utilities.h"
#include "utilities/mesh_utilities.h"
#include "geometries/triangle_3d_3n.h"
#include "utilities/multi_knotspan_boxes_utilities.h"
#include "utilities/integration_points/integration_points_factory.h"
#include "modeler/cube_modeler.h"
//#include "io/io_utilities.h"

//// External includes
#include <fstream>
#include <omp.h>

// TODO: Put enums inside outside test and IntegrationMethod into paramters!
extern std::chrono::duration<double> Elapsed_Time;
/// Functions
void TrIGA::Run(){
  // Obtain discretization of background mesh.
  const double delta_x = std::abs(mParameters.PointA()[0] - mParameters.PointB()[0])/(mParameters.NumberOfElements()[0]);
  const double delta_y = std::abs(mParameters.PointA()[1] - mParameters.PointB()[1])/(mParameters.NumberOfElements()[1]);
  const double delta_z = std::abs(mParameters.PointA()[2] - mParameters.PointB()[2])/(mParameters.NumberOfElements()[2]);

  std::array<double, 3> cube_point_A = {0.0, 0.0, 0.0};
  std::array<double, 3> cube_point_B = {0.0, 0.0, 0.0};
  std::array<double, 3> cube_point_A_param = {0.0, 0.0, 0.0};
  std::array<double, 3> cube_point_B_param = {0.0, 0.0, 0.0};
  const int number_elements_w = mParameters.NumberOfElements()[2];
  const int number_elements_v = mParameters.NumberOfElements()[1];
  const int number_elements_u = mParameters.NumberOfElements()[0];

  const int global_number_of_elements = number_elements_u * number_elements_v * number_elements_w;
  mpElementContainer->reserve(global_number_of_elements);
  int num = 0;
  int trimming_count = 0;
  std::chrono::duration<double> time_intersection{};
  #pragma omp parallel for firstprivate(num, cube_point_A, cube_point_B, cube_point_A_param, cube_point_B_param) schedule(dynamic)
  for( int i = 0; i < global_number_of_elements; ++i){
    // Unroll for loop to enable better parallelization.
    // First walk along rows (x), then columns (y) then into depths (z).
    const int index_in_row_column_plane = i % (number_elements_u*number_elements_v);
    int xx = index_in_row_column_plane % number_elements_u; // row
    int yy = index_in_row_column_plane / number_elements_u; // column
    int zz = i / (number_elements_u*number_elements_v);   // depth

    // Construct bounding box for each element
    cube_point_A[0] = mParameters.PointA()[0] + (xx)*delta_x;
    cube_point_A[1] = mParameters.PointA()[1] + (yy)*delta_y;
    cube_point_A[2] = mParameters.PointA()[2] + (zz)*delta_z;;
    cube_point_B[0] = mParameters.PointA()[0] + (xx+1)*delta_x;
    cube_point_B[1] = mParameters.PointA()[1] + (yy+1)*delta_y;
    cube_point_B[2] = mParameters.PointA()[2] + (zz+1)*delta_z;

    cube_point_A_param = MappingUtilities::FromGlobalToLocalSpace(cube_point_A, mParameters.PointA(), mParameters.PointB());
    cube_point_B_param = MappingUtilities::FromGlobalToLocalSpace(cube_point_B, mParameters.PointA(), mParameters.PointB());

    // Construct element and check status
    std::shared_ptr<Element> tmp_element = std::make_shared<Element>(i+1, cube_point_A_param, cube_point_B_param, mParameters);

    InsideTest::IntersectionStatus status;
    auto cube = CubeModeler::make_cube_3(cube_point_A, cube_point_B);

    if( mEmbeddingFlag ){
      #pragma omp critical // Todo mak this threadsafe!!
      status = mpInsideTest->check_intersection_status_via_element_vertices( *tmp_element);

      // Check if any actually trimmed elements have not been catched previously.
      if( status != InsideTest::Trimmed){

        #pragma omp critical
        trimming_count++;

        if( MeshUtilities::DoIntersect( cube, mPolyhedron) ){
          status = InsideTest::Trimmed;
        }
      }
    }
    else { // If flag is false, consider all knotspans/ elements as inside
      status = InsideTest::Inside;
    }

    bool valid_element = false;

    // Distinguish between trimmed and non-trimmed elements.
    if( status == InsideTest::Trimmed) {


      tmp_element->SetIsTrimmed(true);
      // Create a cubic mesh for each trimmed knotspan
      auto start_time = std::chrono::high_resolution_clock::now();

      auto intersection_mesh = MeshUtilities::GetIntersection(mPolyhedron, cube, mParameters);

      auto end_time = std::chrono::high_resolution_clock::now();

      #pragma omp critical
      time_intersection += end_time - start_time;

      if( intersection_mesh->GetNumberOfCells() > 0 ){
          valid_element = true;
          tmp_element->pSetSurfaceMesh(intersection_mesh);
      }

      if( valid_element ){
        // ExportVolumeMesh(tmp_element->GetSurfaceMesh(), tmp_element->GetId() );
        // Get Gauss-Legendre points in fictitios domain (Only the ones that are outside the material domain)
        // IntegrationPointUtilities::( (*mpInsideTest), tmp_element->GetIntegrationPointsFictitious(),
        //  cube_point_A_param, cube_point_B_param, mParameters.Order()[0], mParameters.Order()[1], mParameters.Order()[2]);

        MomentFitting::ComputeReducedPointsSurfaceIntegral(*tmp_element, mParameters);
        if( tmp_element->GetIntegrationPointsTrimmed().size() == 0 ){
          valid_element = false;
        }
      }

    }
    else if( status == InsideTest::Inside){
      // Get standard gauss legendre points
      if( mParameters.IntegrationMethod() == IntegrationPointFactory::IntegrationMethod::Gauss ){
        IntegrationPointUtilities::CreateGaussLegendrePoints(tmp_element->GetIntegrationPointsInside(),
          cube_point_A_param, cube_point_B_param, mParameters.Order()[0], mParameters.Order()[1], mParameters.Order()[2]);
      }
      else if( mParameters.IntegrationMethod() == IntegrationPointFactory::IntegrationMethod::ReducedGauss1 ){
         IntegrationPointUtilities::CreateGaussLegendrePoints(tmp_element->GetIntegrationPointsInside(),
          cube_point_A_param, cube_point_B_param, mParameters.Order()[0]-1, mParameters.Order()[1]-1, mParameters.Order()[2]-1);
      }
      else if( mParameters.IntegrationMethod() == IntegrationPointFactory::IntegrationMethod::ReducedGauss2 ){
         IntegrationPointUtilities::CreateGaussLegendrePoints(tmp_element->GetIntegrationPointsInside(),
          cube_point_A_param, cube_point_B_param, mParameters.Order()[0]-2, mParameters.Order()[1]-2, mParameters.Order()[2]-2);
      }
      //ExportVolumeMesh(cube, tmp_element->GetId());
      valid_element = true;
    }

    if( valid_element ){
      #pragma omp critical
      mpElementContainer->AddElement(tmp_element); // After this tmp_element is a null_ptr. Is std::moved to container.
    }
  }

  // std::cout << "time_intersection:" << time_intersection.count() << std::endl;
  // std::cout << "trimming_count:  " << trimming_count << std::endl;
  std::cout << "TrIGA :: Elapsed Time: " << Elapsed_Time.count() << std::endl;
  if( (mParameters.IntegrationMethod() != IntegrationPointFactory::IntegrationMethod::Gauss)
        && (mParameters.IntegrationMethod() != IntegrationPointFactory::IntegrationMethod::ReducedGauss1)
          && (mParameters.IntegrationMethod() != IntegrationPointFactory::IntegrationMethod::ReducedGauss2) ){

    //#pragma omp single
    //MultiKnotspanUtilities::ComputeIntegrationPoints(*mpElementContainer, mParameters);

    #pragma omp single
    MultiKnotspanBoxesUtilities::ComputeIntegrationPoints(*mpElementContainer, mParameters);
  }
}
