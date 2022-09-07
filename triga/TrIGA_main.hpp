// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TrIGA_H
#define TrIGA_H

/// CGAL includes
// Domain
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/STL.h>
/// Mesh Processing
#include <CGAL/Polygon_mesh_processing/measure.h>

/// Project includes
#include "io/io_utilities.h"
#include "utilities/intersection_test.h"
#include "geometries/element.h"
#include "utilities/mapping_utilities.h"
#include "utilities/parameters.h"
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
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Mesh_polyhedron_3<K>::type Mesh;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
typedef std::size_t  SizeType;

typedef std::vector<Element> ElementVectorType;

class TrIGA
{
public:

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
            std::cout << "TIBRA :: Start: " << std::endl;

        // Allocate element/knotspans container
        mpElementContainer = std::make_unique<ElementContainer>(mParameters);

        // Read geometry
        if( mEmbeddingFlag ) {
            CGAL::IO::read_STL(mFilename, mPolyhedron);

            // Write Surface Mesh to vtk file if eco_level > 0
            if( mParameters.EchoLevel() > 0){
                IO::WriteMeshToVTK(mPolyhedron, "output/geometry.vtk", true);
            }
            mpIntersectionTest = std::make_unique<IntersectionTest>(mPolyhedron, mParameters.PointA(), mParameters.PointB());
            // Compute volume
            const double volume_global_surface_mesh = CGAL::Polygon_mesh_processing::volume(mPolyhedron);
            if( mParameters.EchoLevel() > 0)
                // std::cout.precision(17);
                std::cout << "Volume of Global Surface Mesh (File: '" << mFilename << "' ): " << volume_global_surface_mesh << std::endl;
        }

        // Start computation
        Run();

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

        if( mParameters.EchoLevel() > 0) {
            // Write vtk files (binary = true)
            IO::WriteElementsToVTK(*mpElementContainer, "output/knotspans.vtk", true);
            IO::WritePointsToVTK(*mpElementContainer, "Trimmed", "output/points_trimmed.vtk", true);
            IO::WritePointsToVTK(*mpElementContainer, "Inside", "output/points_inside.vtk", true);

            std::cout << "TIBRA :: Number of active knotspans: " << mpElementContainer->size() << std::endl;
            std::cout << "TIBRA :: Number of trimmed knotspans: " << number_of_trimmed_elements << std::endl;
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;
        std::cout << "TIBRA :: Elapsed Time: " << elapsed_time.count() << std::endl;
    }

    // Public Member Functions
    std::unique_ptr<ElementContainer> GetElementContainer(){ // Maybe Rename to get ElementContainer
        return std::move(mpElementContainer);
    }

    ElementContainer::ElementVectorPtrType& GetElements(){
        return mpElementContainer->GetElements();
    }

    void ReadWritePostMesh(const std::string& Filename){
        CGAL::IO::read_STL(Filename, mPolyhedronPost);
        IO::WriteMeshToVTK(mPolyhedronPost, "output/results.vtk", true);
    }

    const SurfaceMeshType& GetPostMesh() const {
        return mPolyhedronPost;
    }

private:
    // Private Members
    SurfaceMeshType mPolyhedron;
    SurfaceMeshType mPolyhedronPost;
    std::unique_ptr<IntersectionTest> mpIntersectionTest;
    std::unique_ptr<ElementContainer> mpElementContainer;
    const std::string mFilename;
    const Parameters mParameters;
    const bool mEmbeddingFlag;

    // Private Member Functions
    void Run();
};

#endif // STL_EMBEDDER_H