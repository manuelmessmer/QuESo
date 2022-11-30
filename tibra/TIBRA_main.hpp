// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TIBRA_HPP
#define TIBRA_HPP

/// External includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <chrono>

/// Project includes
#include "io/io_utilities.h"
#include "containers/element.h"
#include "utilities/mapping_utilities.h"
#include "utilities/parameters.h"
#include "containers/element_container.h"
#include "embedding/brep_operator_factory.h"


///@name TIBRA Classes
///@{

////
/**
 * @class  TIBRA
 * @author Manuel Messmer
 * @brief  Main class of TIBRA.
*/
class TIBRA
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t  SizeType;
    typedef std::array<double, 3> PointType;
    typedef std::vector<Element> ElementVectorType;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Constructor. Runs all processes.
    /// @param filename
    /// @param PointA
    /// @param PointB
    /// @param NumberOfElements
    /// @param Order
    /// @param InitialTriangleEdgeLength
    /// @param MinimumNumberOfTriangles
    /// @param MomentFittingResidual
    /// @param PointDistributionFactor
    /// @param IntegrationMethod
    /// @param EchoLevel
    /// @param EmbeddingFlag
    TIBRA(const std::string filename,
                PointType PointA,
                PointType PointB,
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
            // Read mesh
            IO::ReadMeshFromSTL(mTriangleMesh, mFilename.c_str());
            // Write Surface Mesh to vtk file if eco_level > 0
            if( mParameters.EchoLevel() > 0){
                IO::WriteMeshToVTK(mTriangleMesh, "output/geometry.vtk", true);
            }
            // Construct BRepOperator
            mpBRepOperator = BRepOperatorFactory::New(mTriangleMesh);

            // Compute volume
            //const double volume_global_surface_mesh = CGAL::Polygon_mesh_processing::volume(mPolyhedron);
            // if( mParameters.EchoLevel() > 0)
            //     // std::cout.precision(17);
            //     std::cout << "Volume of Global Surface Mesh (File: '" << mFilename << "' ): " << volume_global_surface_mesh << std::endl;
        }

        // Start computation
        Run();

        // Count number of trimmed elements
        SizeType number_of_trimmed_elements = 0;
        std::for_each(mpElementContainer->begin(), mpElementContainer->end(), [&number_of_trimmed_elements] (auto& el_it)
            { if( el_it->IsTrimmed() ) { number_of_trimmed_elements++; } });

        if( mParameters.EchoLevel() > 0) {
            // Write vtk files (binary = true)
            IO::WriteElementsToVTK(*mpElementContainer, "output/knotspans.vtk", true);
            IO::WritePointsToVTK(*mpElementContainer, "Trimmed", "output/points_trimmed.vtk", true);
            IO::WritePointsToVTK(*mpElementContainer, "Inside", "output/points_inside.vtk", true);

            std::cout << "TIBRA :: Number of active knotspans: " << mpElementContainer->size() << std::endl;
            std::cout << "TIBRA :: Number of trimmed knotspans: " << number_of_trimmed_elements << std::endl;

            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_time = end_time - start_time;
            std::cout << "TIBRA :: Elapsed Time: " << elapsed_time.count() << std::endl;
        }
    }

    /// Copy Constructor
    TIBRA(const TIBRA &m) = delete;
    /// Copy Assignement
    TIBRA & operator= (const TIBRA &) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Get all active elements.
    /// @return const Reference to ElementVectorPtrType
    const ElementContainer::ElementVectorPtrType& GetElements() const {
        return mpElementContainer->GetElements();
    }

    /// @brief Reads Filename and writes mesh to output/results.vtk
    /// @param Filename
    void ReadWritePostMesh(const std::string& rFilename) {
        IO::ReadMeshFromSTL(mTriangleMeshPost, rFilename.c_str());
        IO::WriteMeshToVTK(mTriangleMeshPost, "output/results.vtk", true);
    }

    /// @brief  Get mesh for prosptrocessing
    /// @return const Reference to TriangleMesh
    const TriangleMesh& GetPostMesh() const {
        return mTriangleMeshPost;
    }
    ///@}

private:

    ///@name Private Members Variables
    ///@{
    TriangleMesh mTriangleMesh;
    TriangleMesh mTriangleMeshPost;
    std::unique_ptr<BRepOperatorBase> mpBRepOperator;
    std::unique_ptr<ElementContainer> mpElementContainer;
    const std::string mFilename;
    const Parameters mParameters;
    const bool mEmbeddingFlag;
    ///@}

    ///@name Private Member Operations
    ///@{

    /// @brief Run TIBRA
    void Run();
    ///@}
};

///@}

#endif // TIBRA_HPP