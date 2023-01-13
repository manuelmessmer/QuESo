// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TIBRA_INCLUDE_H
#define TIBRA_INCLUDE_H

/// STL includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>

/// Project includes
#include "containers/element_container.hpp"
#include "embedding/brep_operator_factory.h"
#include "containers/element.hpp"
#include "io/io_utilities.h"
#include "utilities/mapping_utilities.h"
#include "utilities/mesh_utilities.h"
#include "utilities/parameters.h"
#include "utilities/timer.hpp"
#include "utilities/logger.hpp"

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
        Timer timer{};
        TIBRA_INFO_IF(mParameters.EchoLevel() > 0) << "\nTIBRA ------------------------------------------ START" << std::endl;

        // Allocate element/knotspans container
        mpElementContainer = MakeUnique<ElementContainer>(mParameters);

        // Read geometry
        double volume_brep = 0.0;
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
            volume_brep = MeshUtilities::VolumeOMP(mTriangleMesh);

            TIBRA_INFO_IF(mParameters.EchoLevel() > 0) << "Read file: '" << r_filename << "'\n";
            TIBRA_INFO_IF(mParameters.EchoLevel() > 0) << "Volume of B-Rep model: " << volume_brep << '\n';
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

            TIBRA_INFO << "Number of active knotspans: " << mpElementContainer->size() << std::endl;
            TIBRA_INFO << "Number of trimmed knotspans: " << number_of_trimmed_elements << std::endl;

            if( mParameters.EchoLevel() > 1 ) {
                const double volume_ips = mpElementContainer->GetVolumeOfAllIPs();
                TIBRA_INFO << "The computed quadrature represents " << volume_ips/volume_brep * 100
                    << "% of the volume of the BRep model.\n";
            }

            TIBRA_INFO << "Elapsed time: " << timer.Measure() << std::endl;
            TIBRA_INFO << "TIBRA ------------------------------------------- END\n" << std::endl;
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

    ///@}
    ///@name Temporary operations to perform PosProcessing after Kratos Analysis
    ///      This will be moved to TriangleMesh.
    ///@{

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
    Unique<BRepOperatorBase> mpBRepOperator;
    Unique<ElementContainer> mpElementContainer;
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

#endif // TIBRA_INCLUDE_H