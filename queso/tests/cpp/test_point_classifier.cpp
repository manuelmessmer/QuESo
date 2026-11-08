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

//// External includes
#include <doctest/doctest.h>
//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/mesh_operator.h"
#include "queso/io/io_utilities.h"

#include "queso/tests/cpp/helpers/test_data.hpp"

namespace queso::testing {

    TEST_SUITE("PointClassifierTestSuite")
    {

        TEST_CASE("PointClassifierTestSuite - Cylinder Point Classifier Test")
        {


            TriangleMesh triangle_mesh{};
            // Read mesh from STL file
            std::string base_dir = TestDataDirectory();
            IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/primitives/cylinder.stl");

            std::vector<PointType> rPoints{};
            rPoints.reserve(167620);
            for (double x = -1.5; x <= 1.5; x += 0.09) {
                for (double y = -1.5; y <= 1.5; y += 0.09) {
                    for (double z = -1; z <= 12; z += 0.09) { rPoints.push_back({ x, y, z }); }
                }
            }

            // Instantiate brep_operator
            embedding::MeshOperator classifier(
                triangle_mesh.View(), GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 2.0 })
            );

            std::vector<bool> result(rPoints.size(), false);
            int count = 0;
            for (auto& point : rPoints) {
                if (classifier.IsInside(point)) { result[count] = true; }
                count++;
            }

            for (IndexType i = 0; i < result.size(); ++i) {
                double radius = std::sqrt(rPoints[i][0] * rPoints[i][0] + rPoints[i][1] * rPoints[i][1]);
                if (radius < 1.0 && rPoints[i][2] > 0.0 && rPoints[i][2] < 10.0) {
                    CHECK((result)[i]);
                } else {
                    CHECK_FALSE((result)[i]);
                }
            }
        }

        TEST_CASE("PointClassifierTestSuite - Cube Point Classifier Test")
        {


            TriangleMesh triangle_mesh{};
            // Read mesh from STL file
            std::string base_dir = TestDataDirectory();
            IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/primitives/cube_with_cavity.stl");

            std::vector<PointType> rPoints{};
            // rPoints.reserve(117900);
            for (double x = -1.5; x <= 1.5 + 1e-15; x += 0.15) {
                for (double y = -1.5; y <= 1.5 + 1e-15; y += 0.15) {
                    for (double z = -1.5; z <= 1.5 + 1e-15; z += 0.15) { rPoints.push_back({ x, y, z }); }
                }
            }

            // Instantiate brep_operator
            embedding::MeshOperator classifier(
                triangle_mesh.View(), GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 2.0 })
            );

            std::vector<bool> result(rPoints.size(), false);
            int count = 0;
            for (auto& point : rPoints) {
                if (classifier.IsInside(point)) { result[count] = true; }
                count++;
            }

            for (IndexType i = 0; i < result.size(); ++i) {
                double radius = std::sqrt(
                    rPoints[i][0] * rPoints[i][0] + rPoints[i][1] * rPoints[i][1] + rPoints[i][2] * rPoints[i][2]
                );
                if (radius <= 1.0) {
                    CHECK_FALSE((result)[i]);
                } else {
                    if (rPoints[i][0] <= -1.5 + 1e-14 || rPoints[i][0] >= 1.5 - 1e-14 || rPoints[i][1] <= -1.5 + 1e-14
                        || rPoints[i][1] >= 1.5 - 1e-14 || rPoints[i][2] <= -1.5 + 1e-14
                        || rPoints[i][2] >= 1.5 - 1e-14) {

                        CHECK_FALSE((result)[i]);
                    } else {
                        CHECK((result)[i]);
                    }
                }
            }

            // Check if point on boundary is not inside.
            // Note the stl has a mesh at (1.0, 0.0, 0.0)
            CHECK_FALSE(classifier.IsInside(PointType{ 1.0, 0.0, 0.0 }));
        }

        TEST_CASE("PointClassifierTestSuite - Elephant Point Classifier Test")
        {


            TriangleMesh triangle_mesh{};
            // Read mesh from STL file
            std::string base_dir = TestDataDirectory();
            IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/primitives/elephant.stl");

            // Instantiate brep_operator
            embedding::MeshOperator classifier(
                triangle_mesh.View(), GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 })
            );

            std::vector<PointType> rPoints{};
            rPoints.reserve(84000);
            for (double x = -0.4; x <= 0.4; x += 0.02) {
                for (double y = -0.6; y <= 0.6; y += 0.02) {
                    for (double z = -0.35; z <= 0.35; z += 0.02) { rPoints.push_back({ x, y, z }); }
                }
            }

            std::vector<bool> result(rPoints.size(), false);
            int count = 0;
            for (auto& point : rPoints) {
                if (classifier.IsInside(point)) { result[count] = true; }
                count++;
            }

            std::vector<bool> result_ref{};
            // Read reference results from file
            std::string line;
            std::ifstream myfile(base_dir + "/embedding/inside_outside_elephant.txt");
            if (myfile.is_open()) {
                while (getline(myfile, line)) { result_ref.push_back(std::stoi(line)); }
                myfile.close();
            }

            // std::ofstream myfile;
            // myfile.open ("test.txt");
            // Compare results
            for (IndexType i = 0; i < result.size(); ++i) { CHECK((result[i]) == (result_ref[i])); }
            // myfile.close();
        }


        TEST_CASE("PointClassifierTestSuite - Bunny Point Classifier Test")
        {


            TriangleMesh triangle_mesh{};
            // Read mesh from STL file
            std::string base_dir = TestDataDirectory();
            IO::ReadMeshFromSTL(triangle_mesh, base_dir + "/primitives/stanford_bunny.stl");

            // Instantiate brep_operator
            embedding::MeshOperator classifier(
                triangle_mesh.View(), GeometryTolerance::FromScale({ .length_scale = 10.0, .coordinate_scale = 115.0 })
            );

            std::vector<PointType> rPoints{};
            rPoints.reserve(138600);
            for (double x = -24; x <= 85; x += 2) {
                for (double y = -43; y <= 46; y += 2) {
                    for (double z = 5; z <= 115; z += 2) { rPoints.push_back({ x, y, z }); }
                }
            }

            std::vector<bool> result(rPoints.size(), false);
            int count = 0;
            for (auto& point : rPoints) {
                if (classifier.IsInside(point)) { result[count] = true; }
                count++;
            }

            std::vector<bool> result_ref{};
            // Read reference results from file
            std::string line;
            std::ifstream myfile(base_dir + "/embedding/inside_outside_bunny.txt");
            if (myfile.is_open()) {
                while (getline(myfile, line)) { result_ref.push_back(std::stoi(line)); }
                myfile.close();
            }

            // std::ofstream myfile;
            // myfile.open ("test.txt");
            // Compare results
            for (IndexType i = 0; i < result.size(); ++i) { CHECK((result[i]) == (result_ref[i])); }
            // myfile.close();
        }

    }  // TEST_SUITE

} // namespace queso::testing
