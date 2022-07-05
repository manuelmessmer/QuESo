// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MOMENT_FITTING_UTILITIES_INCLUDE_H
#define MOMENT_FITTING_UTILITIES_INCLUDE_H

// External includes
#include <vector>
#include <array>

// Project includes
#include "geometries/element.h"
#include "utilities/parameters.h"

typedef Element::IntegrationPointVectorType IntegrationPointVectorType;

namespace MomentFitting{

double ComputeReducedPointsSurfaceIntegral(Element& rElement, const int PointDistributionFactor, const Parameters& rParam);

void ComputeReducedPointsSurfaceIntegral(Element& rElement, const Parameters& rParam);

void DistributeIntegrationPoints(Element& rElement, IntegrationPointVectorType& rIntegrationPoint, const int PointDistributionFactor, const Parameters& rParam);

double f_x(double x, int order);

double f_x_integral(double x, int order);

double p_n(double x, int order);

double f_x(double x, int order, double a, double b);

double f_x_integral(double x, int order, double a, double b);

} // End Namespace

#endif // MOMENT_FITTING_UTILITIES_INCLUDE_H