//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

//// Project includes
#include "queso/quadrature/moment_fitting_assembly.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/integration_point.hpp"

namespace queso::quadrature::moment_fitting::detail {

// Explicit template instantiations for QuESo's native integration-point vector types.
template std::vector<double> ComputeConstantTerms<std::vector<BoundaryIntegrationPoint>>(
    const std::vector<BoundaryIntegrationPoint>& rBoundaryIPs,
    const BoundingBoxType& rBoundsGlobal,
    const IntegrationOrderInfo& rOrderInfo
);

template void AssembleFittingMatrix<std::vector<IntegrationPoint>>(
    NNLS::MatrixType& rFittingMatrix,
    const std::vector<IntegrationPoint>& rIntegrationPoints,
    const BoundingBoxType& rBoundsParam,
    const IntegrationOrderInfo& rOrderInfo
);

}  // namespace queso::quadrature::moment_fitting::detail
