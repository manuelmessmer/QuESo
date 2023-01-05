// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TIBRA_HPP
#define TIBRA_HPP

/// STL includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <chrono>

/// Project includes
#include "containers/element_container.h"
#include "embedding/brep_operator_factory.h"
#include "containers/element.h"
#include "io/io_utilities.h"
#include "utilities/mapping_utilities.h"
#include "utilities/mesh_utilities.h"
#include "utilities/parameters.h"

namespace tibra {

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

    typedef std::vector<Element> ElementVectorType;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Constructor. Runs all processes.
    TIBRA(const Parameters& rParameters ) : mParameters(rParameters)
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        if( mParameters.EchoLevel() > 0)
            std::cout << "TIBRA :: Start: " << std::endl;

        // Allocate element/knotspans container
        mpElementContainer = std::make_unique<ElementContainer>(mParameters);

        // Read geometry
        if( mParameters.Get<bool>("embedding_flag") ) {
            // Read mesh
            const auto& r_filename = mParameters.Get<std::string>("input_filename");
            IO::ReadMeshFromSTL(mTriangleMesh, r_filename.c_str());
            // Write Surface Mesh to vtk file if eco_level > 0
            if( mParameters.EchoLevel() > 0){
                IO::WriteMeshToVTK(mTriangleMesh, "output/geometry.vtk", true);
            }
            // Construct BRepOperator
            mpBRepOperator = BRepOperatorFactory::New(mTriangleMesh, mParameters);

            // Compute volume
            const double volume = MeshUtilities::Volume(mTriangleMesh);
            if( mParameters.EchoLevel() > 0){
                std::cout << "TIBRA :: Read file: '" << r_filename << "'\n";
                std::cout << "TIBRA :: Volume of input model: " << volume << '\n';
            }
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
            IO::WritePointsToVTK(*mpElementContainer, "All", "output/integration_points.vtk", true);

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
    void ReadWritePostMesh() {
        const auto& r_filename = mParameters.Get<std::string>("postprocess_filename");
        IO::ReadMeshFromSTL(mTriangleMeshPost, r_filename.c_str());
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
    const Parameters mParameters;
    ///@}

    ///@name Private Member Operations
    ///@{

    /// @brief Run TIBRA
    void Run();
    ///@}
};

///@}

} // End namespace tibra

#endif // TIBRA_HPP