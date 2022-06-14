// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TrIGA_H
#define TrIGA_H

// VTK includes
#include <vtkSTLReader.h>
#include <vtkMassProperties.h>
#include <vtkPolyDataNormals.h>

/// Project includes
//#include "io/vtk_utilities.h"
#include "utilities/inside_test.h"
#include "geometries/element.h"
#include "utilities/mapping_utilities.h"
#include "utilities/parameters.h"
#include "io/io_utilities.h"
#include "geometries/element_container.h"

/// External includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <chrono>

/// Type Definitions

class TrIGA
{
public:

    // Typedefs
    typedef std::size_t  SizeType;
    typedef std::size_t  IndexType;
    typedef std::vector<Element> ElementVectorType;

    // Constructor
    TrIGA(const std::string filename,
                std::array<double, 3> PointA,
                std::array<double, 3> PointB,
                std::array<int, 3> NumberOfElements,
                std::array<int, 3> Order,
                double InitialTriangleEdgeLength,
                double MinimumNumberOfTriangles,
                double MomentFittingResidual,
                double PointDistributionFactor,
                std::string IntegrationMethod,
                int EchoLevel,
                bool EmbeddingFlag = true) :
            mFilename(filename),
            mEmbeddingFlag(EmbeddingFlag),
            mParameters(PointA, PointB, NumberOfElements, Order,
                        InitialTriangleEdgeLength, MinimumNumberOfTriangles, MomentFittingResidual,
                        PointDistributionFactor, IntegrationMethod, EchoLevel)
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        if( mParameters.EchoLevel() > 0)
            std::cout << "TrIGA :: Start: " << std::endl;

        // Allocate element/knotspans container
        mpElementContainer = std::make_unique<ElementContainer>(mParameters);

        // Read geometry
        if( mEmbeddingFlag ) {
            // Read STL
            vtkSmartPointer<vtkSTLReader> sr = vtkSmartPointer<vtkSTLReader>::New();
            sr->SetFileName(filename.c_str());
            sr->Update();
            mPolyhedron = sr->GetOutput();

            vtkSmartPointer<vtkPolyDataNormals> dataset = vtkSmartPointer<vtkPolyDataNormals>::New();
            dataset->SetInputData(mPolyhedron);
            dataset->ComputePointNormalsOff();
            dataset->ComputeCellNormalsOn();
            dataset->Update();

            mPolyhedron->ShallowCopy( dataset->GetOutput() );

            if (!mPolyhedron)
            {
                throw  std::runtime_error("Mesh is not a valid stl file. File: " + mFilename + " has been tried to read. \n");
            }

            // Write Surface Mesh to vtk file if echo_level > 0
            if( mParameters.EchoLevel() > 0){
                IO_Utilities::WriteVTK(mPolyhedron, "output/geometry.vtk");
            }
            mpInsideTest = std::make_unique<InsideTest>(mPolyhedron, mParameters.PointA(), mParameters.PointB());

            // Compute Volume if echo_level > 0
            if( mParameters.EchoLevel() > 0){
                vtkSmartPointer<vtkMassProperties> mass_properties = vtkSmartPointer<vtkMassProperties>::New();
                mass_properties->SetInputData(mPolyhedron);
                mass_properties->Update();
                //std::cout.precision(17);
                std::cout << "Volume of Global Surface Mesh (File: '" << mFilename << "' ): " << mass_properties->GetVolume() << std::endl;
            }
        }

        // Start computation
        this->Run();

        // Compute number of trimmed elements
        auto element_it_begin = mpElementContainer->begin();
        SizeType number_of_trimmed_elements = 0;
        for( int i = 0; i < mpElementContainer->size(); ++i ){
            auto element_it = element_it_begin + i;
            if( (*element_it)->IsTrimmed() ) {
                number_of_trimmed_elements++;

                const int number_points = (*element_it)->GetIntegrationPointsTrimmed().size();
                //TODO: Make this generic!
                // if( number_points < 8 || number_points > 27 ) { //This is valid for p=2.
                //     std::stringstream error_message;
                //     error_message << "Inappropriate number of integration points in trimmed element with ID: ";
                //     error_message << (*element_it)->GetId() << ". Number of Integration Points: " << number_points << ".\n";
                //     error_message << "Targeted number of integration points: 8 < x < 27. Note: These values are hardcoded for ansatz order p=2." << std::endl;
                //     throw std::runtime_error(error_message.str());
                // }
            }
            else {
                const int number_points = (*element_it)->GetIntegrationPointsInside().size();
                // if( number_points < 8 || number_points > 27){ // This is valid for p=2.
                //     std::stringstream error_message;
                //     error_message << "Inappropriate number of integration points in non-trimmed element with ID: ";
                //     error_message << (*element_it)->GetId() << ". Number of Integration Points: " << number_points << ".\n";
                //     error_message << "Targeted number of integration points: 8 < x < 27. Note: These values are hardcoded for ansatz order p=2." << std::endl;
                //     throw std::runtime_error(error_message.str());
                // }
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;
        if( mParameters.EchoLevel() > 0) {
            std::cout << "TrIGA :: Number of active knotspans: " << mpElementContainer->size() << std::endl;
            std::cout << "TrIGA :: Number of trimmed knotspans: " << number_of_trimmed_elements << std::endl;
            std::cout << "TrIGA :: Elapsed Time: " << elapsed_time.count() << std::endl;
        }
    }

    // Public Member Functions
    std::unique_ptr<ElementContainer> GetElements(){ // Maybe Rename to get ElementContainer
        return std::move(mpElementContainer);
    }

    ElementContainer::ElementVectorPtrType& ExportElements(){
        return mpElementContainer->GetElements();
    }

private:
    // Private Members
    vtkSmartPointer<vtkPolyData> mPolyhedron;
    std::unique_ptr<InsideTest> mpInsideTest;
    std::vector<Element::IntegrationPointPtrType> mPoints{};
    std::unique_ptr<ElementContainer> mpElementContainer;
    const std::string mFilename;
    const Parameters mParameters;
    const bool mEmbeddingFlag;

    // Private Member Functions
    void Run();

    // std::vector<IntegrationPoint*> GetPoints(){
    //     const auto element_itr_begin = mpElementContainer->begin();
    //     for( int i = 0; i < mpElementContainer->size(); ++i){
    //         auto element_itr = *(element_itr_begin + i);

    //         mPoints.insert( mPoints.end(), element_itr->GetIntegrationPointsTrimmed().begin(), element_itr->GetIntegrationPointsTrimmed().end() );
    //     }
    // }
};

#endif // TrIGA_H