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

#pragma once

//// STL includes
#include <cmath>

//// Project includes
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace queso::quadrature::tensor_product {

/// @brief Parameters for tensor-product quadrature rules on a single non-trimmed element.
/// @details Available quadrature rules: {Gauss, Gauss_Reduced1, Gauss_Reduced2}.
struct Parameters
{
    Vector3i integration_order;
    IntegrationMethodType method = IntegrationMethod::gauss;
};

namespace detail {

    template<typename TIntegrationPointVectorType>
    void ComputeImpl(
        TIntegrationPointVectorType& rIntegrationPoints,
        const BoundingBoxType& rBoundsParam,
        const Parameters& rParameters
    )
    {
        const auto& r_ip_list_u =
            IntegrationPointFactory1D::GetGauss(rParameters.integration_order[0], rParameters.method);
        const auto& r_ip_list_v =
            IntegrationPointFactory1D::GetGauss(rParameters.integration_order[1], rParameters.method);
        const auto& r_ip_list_w =
            IntegrationPointFactory1D::GetGauss(rParameters.integration_order[2], rParameters.method);

        const SizeType n_point = r_ip_list_u.size() * r_ip_list_v.size() * r_ip_list_w.size();
        rIntegrationPoints.clear();
        rIntegrationPoints.reserve(n_point);

        const auto& lower = rBoundsParam.lower;
        const auto& upper = rBoundsParam.upper;

        const double lower_u = lower[0];
        const double lower_v = lower[1];
        const double lower_w = lower[2];

        const double length_u = std::abs(upper[0] - lower_u);
        const double length_v = std::abs(upper[1] - lower_v);
        const double length_w = std::abs(upper[2] - lower_w);

        for (const auto& r_ip_u : r_ip_list_u) {
            const double u = lower_u + length_u * r_ip_u[0];
            const double weight_u = r_ip_u[1] * length_u;

            for (const auto& r_ip_v : r_ip_list_v) {
                const double v = lower_v + length_v * r_ip_v[0];
                const double weight_uv = weight_u * r_ip_v[1] * length_v;

                for (const auto& r_ip_w : r_ip_list_w) {
                    const double w = lower_w + length_w * r_ip_w[0];
                    rIntegrationPoints.emplace_back(u, v, w, weight_uv * r_ip_w[1] * length_w);
                }
            }
        }
    }
}  // namespace detail

/// @brief Assembles tensor-product quadrature points on the given parametric bounds.
/// @note This function clears rIntegrationPoints before adding the new ones.
/// @param[out] rIntegrationPoints Integration point container.
/// @param rBoundsParam Bounds in parametric space.
/// @param rParameters Quadrature parameters.
template<typename TIntegrationPointVectorType>
void Compute(
    TIntegrationPointVectorType& rIntegrationPoints,
    const BoundingBoxType& rBoundsParam,
    const Parameters& rParameters
)
{ detail::ComputeImpl(rIntegrationPoints, rBoundsParam, rParameters); }

/// @brief Assembles tensor-product quadrature points on an element.
/// @note This function clears the element's existing integration points before adding the new ones.
/// @param rElement Element whose parametric cell bounds are used.
/// @param rParameters Quadrature parameters.
template<typename TElementType>
void Compute(TElementType& rElement, const Parameters& rParameters)
{
    Compute(
        rElement.GetIntegrationPoints(), rElement.template GetCellBounds<CoordinateSpace::parametric>(), rParameters
    );
}

}  // namespace queso::quadrature::tensor_product
