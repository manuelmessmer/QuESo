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

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/io/io_utilities.h"
#include "queso/embedding/brep_operator.h"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( ElementClassifierTestSuite )

BOOST_AUTO_TEST_CASE(TouchElementCubeTest) {
    QuESo_INFO << "Testing :: Test Classify Elements :: Touch Cube" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cube_with_cavity.stl");

    // Instatiate brep_operator
    BRepOperator brep_operator(triangle_mesh);

    Vector3d lower_bound = {-2, -2, -2};
    Vector3d upper_bound = {-1.5, 2, 2};
    // Touch from outside with tolerance=0.0 is trimmed.
    QuESo_CHECK_EQUAL( brep_operator.GetIntersectionState(lower_bound, upper_bound, 0.0), IntersectionState::trimmed );
    // Touch from outside with tolerance>0.0 is outside.
    QuESo_CHECK_EQUAL( brep_operator.GetIntersectionState(lower_bound, upper_bound, 1e-8), IntersectionState::outside );

    lower_bound = {-1.5, -1.5, -1.5};
    upper_bound = {-1.4, -1.4, -1.4};
    // Touch from inside with tolerance=0.0 is trimmed.
    QuESo_CHECK_EQUAL( brep_operator.GetIntersectionState(lower_bound, upper_bound, 0.0), IntersectionState::trimmed );
    // Touch from inside with tolerance>0.0 is inside.
    QuESo_CHECK_EQUAL( brep_operator.GetIntersectionState(lower_bound, upper_bound, 1e-8), IntersectionState::inside );
}

BOOST_AUTO_TEST_CASE(CylinderElementClassifierTest) {
    QuESo_INFO << "Testing :: Test Classify Elements :: Cylinder" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cylinder.stl");

    Settings settings;
    auto& r_grid_settings = settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-1.5, -1.5, -1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{1.5, 1.5, 12.0});
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{-1.0, -1-0, 1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 1.0, 1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{30, 30, 130});

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh);
    GridIndexer grid_indexer(settings);

    const SizeType num_of_elements = grid_indexer.NumberOfElements();
    std::vector<IndexType> result{};
    result.reserve(num_of_elements);
    std::vector<std::pair<PointType, PointType>> boxes{};
    for(IndexType i = 0; i < num_of_elements; ++i) {
        auto bounding_box = grid_indexer.GetBoundingBoxXYZFromIndex(i);
        result.push_back( static_cast<IndexType>(brep_operator.GetIntersectionState(bounding_box.first, bounding_box.second)) );

    }

    // Get flood fill solution
    const auto p_classification = brep_operator.pGetElementClassifications(settings);

    QuESo_CHECK_EQUAL( result.size(), num_of_elements);
    std::ifstream myfile("queso/tests/cpp_tests/results/element_classifier_cylinder.txt");
    std::string line;
    for( IndexType i = 0; i < result.size(); ++i){
        getline (myfile, line);
        QuESo_CHECK_EQUAL( result[i], static_cast<IndexType>(std::stoi(line)) );
        // Now test against flood fill solution.
        QuESo_CHECK_EQUAL( result[i], static_cast<IndexType>((*p_classification)[i]) );
    }
    myfile.close();
}

BOOST_AUTO_TEST_CASE(CubeElementClassifierTest) {
    QuESo_INFO << "Testing :: Test Classify Elements :: Cube with cavity" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cube_with_cavity.stl");

    Settings settings;
    auto& r_grid_settings = settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-1.5, -1.5, -1.5});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{1.5, 1.5, 1.5});
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{-1.0, -1-0, 1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 1.0, 1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{20, 20, 20});

    // Instatiate brep_operator
    BRepOperator brep_operator(triangle_mesh);
    GridIndexer grid_indexer(settings);

    const SizeType num_of_elements = grid_indexer.NumberOfElements();
    std::vector<IndexType> result{};
    result.reserve(num_of_elements);
    for( IndexType i = 0; i < num_of_elements; ++i ){
        const auto bounding_box = grid_indexer.GetBoundingBoxXYZFromIndex(i);
        result.push_back( static_cast<IndexType>(brep_operator.GetIntersectionState(bounding_box.first, bounding_box.second)) );
    }

    // Get flood fill solution
    const auto p_classification = brep_operator.pGetElementClassifications(settings);

    QuESo_CHECK_EQUAL( result.size(), num_of_elements);
    std::ifstream myfile("queso/tests/cpp_tests/results/element_classifier_cube.txt");
    std::string line;
    for( IndexType i = 0; i < result.size(); ++i){
        getline (myfile, line);
        QuESo_CHECK_EQUAL( result[i], static_cast<IndexType>(std::stoi(line)) );
        // Now test against flood fill solution.
        QuESo_CHECK_EQUAL( result[i], static_cast<IndexType>((*p_classification)[i]) );
    }
    myfile.close();
}

BOOST_AUTO_TEST_CASE(ElephantElementClassifierTest) {
    QuESo_INFO << "Testing :: Test Classify Elements :: Elephant" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/elephant.stl");

    Settings settings;
    auto& r_grid_settings = settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-0.4, -0.6, -0.35});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{0.4, 0.6, 0.35});
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{-1.0, -1-0, 1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 1.0, 1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{16, 24, 14});

    // Instatiate brep_operator
    BRepOperator brep_operator(triangle_mesh);
    GridIndexer grid_indexer(settings);

    std::vector<IndexType> result{};
    result.reserve(5376);
    for( IndexType i = 0; i < 5376; ++i ){
        const auto bounding_box = grid_indexer.GetBoundingBoxXYZFromIndex(i);
        result.push_back( static_cast<IndexType>(brep_operator.GetIntersectionState(bounding_box.first, bounding_box.second)) );
    }

    // Get flood fill solution
    const auto p_classification = brep_operator.pGetElementClassifications(settings);

    QuESo_CHECK_EQUAL( result.size(), 5376);
    std::ifstream myfile("queso/tests/cpp_tests/results/element_classifier_elephant.txt");
    std::string line;
    for( IndexType i = 0; i < result.size(); ++i){
        getline (myfile, line);
        QuESo_CHECK_EQUAL( result[i], static_cast<IndexType>( std::stoi(line)) );
        // Now test against flood fill solution.
        QuESo_CHECK_EQUAL( result[i], static_cast<IndexType>((*p_classification)[i]) );
    }
    myfile.close();
}

BOOST_AUTO_TEST_CASE(BunnyElementClassifierTest) {
    QuESo_INFO << "Testing :: Test Classify Elements :: Bunny" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/stanford_bunny.stl");

    Settings settings;
    auto& r_grid_settings = settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-24, -43, 5});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{85, 46, 115});
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{-1.0, -1-0, 1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 1.0, 1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{36, 30, 40});

    // Instatiate brep_operator
    BRepOperator brep_operator(triangle_mesh);
    GridIndexer grid_indexer(settings);

    const SizeType num_of_elements = grid_indexer.NumberOfElements();
    std::vector<IndexType> result{};
    result.reserve(num_of_elements);
    for( IndexType i = 0; i < num_of_elements; ++i){
        const auto bounding_box = grid_indexer.GetBoundingBoxXYZFromIndex(i);
        result.push_back( static_cast<IndexType>(brep_operator.GetIntersectionState(bounding_box.first, bounding_box.second )) );
    }

    // Get flood fill solution
    const auto p_classification = brep_operator.pGetElementClassifications(settings);

    QuESo_CHECK_EQUAL( result.size(), num_of_elements);
    std::ifstream myfile("queso/tests/cpp_tests/results/element_classifier_bunny.txt");
    std::string line;
    for( IndexType i = 0; i < result.size(); ++i){
        getline (myfile, line);
        QuESo_CHECK_EQUAL( result[i], static_cast<IndexType>(std::stoi(line)) );
        // Now test against flood fill solution.
        QuESo_CHECK_EQUAL( result[i], static_cast<IndexType>((*p_classification)[i]) );
    }
    myfile.close();
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso