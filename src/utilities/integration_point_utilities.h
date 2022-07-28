// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef INTEGRATION_POINT_UTILITIES_H
#define INTEGRATION_POINT_UTILITIES_H

#include <vector>
#include <array>

// Project includes
#include "geometries/integration_point.h"
#include "utilities/inside_test.h"

namespace IntegrationPointUtilities {

typedef std::vector<IntegrationPoint> IntegrationPointType;
typedef size_t SizeType;

void IntegrationPoints3D(
        InsideTest& rInsideTest,
        IntegrationPointType& rIntegrationPoints,
        SizeType PointsInU, SizeType PointsInV, SizeType PointsInW,
        double U0, double U1, double V0, double V1, double W0, double W1);

void IntegrationPoints3D(
        IntegrationPointType& rIntegrationPoints,
        SizeType PointsInU, SizeType PointsInV, SizeType PointsInW,
        double U0, double U1, double V0, double V1, double W0, double W1);

void CreateGaussLegendrePoints(
        InsideTest& rInsideTest,
        IntegrationPointType& rIntegrationPoints,
        std::array<double,3> point_A,
        std::array<double,3> point_B,
        SizeType OrderU,
        SizeType OrderV,
        SizeType OrderW);

void CreateGaussLegendrePoints(
        IntegrationPointType& rIntegrationPoints,
        std::array<double,3> point_A,
        std::array<double,3> point_B,
        SizeType OrderU,
        SizeType OrderV,
        SizeType OrderW);

} // End IntegrationPointUtilities

#endif // INTEGRATION_POINT_UTILITIES_H