// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "containers/triangle_mesh.hpp"
#include "io/io_utilities.h"
#include "embedding/brep_operator.h"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( ElementClassifierTestSuite )

BOOST_AUTO_TEST_CASE(TouchElementCubeTest) {
    TIBRA_INFO << "Testing :: Test Classify Elements :: Touch Cube" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

    Parameters params( {Component("lower_bound_xyz", PointType(0.0, 0.0, 0.0)),
                        Component("upper_bound_xyz", PointType(1.0, 1.0, 1.0)),
                        Component("lower_bound_uvw", PointType(0.0, 0.0, 0.0)),
                        Component("upper_bound_uvw", PointType(1.0, 1.0, 1.0)),
                        Component("number_of_elements", Vector3i(1, 1, 1)) });

    // Instatiate brep_operator
    BRepOperator brep_operator(triangle_mesh, params);

    Vector3d lower_bound = {-2, -2, -2};
    Vector3d upper_bound = {-1.5, 2, 2};
    // Touch from outside with tolerance=0.0 is trimmed.
    BOOST_CHECK_EQUAL( brep_operator.GetIntersectionState(lower_bound, upper_bound, 0.0), IntersectionStatus::Trimmed );
    // Touch from outside with tolerance>0.0 is outside.
    BOOST_CHECK_EQUAL( brep_operator.GetIntersectionState(lower_bound, upper_bound, 1e-8), IntersectionStatus::Outside );

    lower_bound = {-1.5, -1.5, -1.5};
    upper_bound = {-1.4, -1.4, -1.4};
    // Touch from inside with tolerance=0.0 is trimmed.
    BOOST_CHECK_EQUAL( brep_operator.GetIntersectionState(lower_bound, upper_bound, 0.0), IntersectionStatus::Trimmed );
    // Touch from inside with tolerance>0.0 is inside.
    BOOST_CHECK_EQUAL( brep_operator.GetIntersectionState(lower_bound, upper_bound, 1e-8), IntersectionStatus::Inside );
}

BOOST_AUTO_TEST_CASE(CylinderElementClassifierTest) {
    TIBRA_INFO << "Testing :: Test Classify Elements :: Cylinder" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    Parameters params( {Component("lower_bound_xyz", PointType(-1.5, -1.5, -1.0)),
                        Component("upper_bound_xyz", PointType(1.5, 1.5, 12.0)),
                        Component("lower_bound_uvw", PointType(-1.5, -1.5, -1.0)),
                        Component("upper_bound_uvw", PointType(1.5, 1.5, 12.0)),
                        Component("number_of_elements", Vector3i(30, 30, 130)) });

    // Instantiate brep_operator
    BRepOperator brep_operator(triangle_mesh, params);
    Mapper mapper(params);

    std::vector<IndexType> result{};
    result.reserve(117000);
    std::vector<std::pair<PointType, PointType>> boxes{};
    for(IndexType i = 0; i < 117000; ++i) {
        auto bounding_box = mapper.GetBoundingBoxXYZFromIndex(i);
        result.push_back( brep_operator.GetIntersectionState(bounding_box.first, bounding_box.second) );

    }

    // Get flood fill solution
    const auto p_classification = brep_operator.pGetElementClassifications();

    BOOST_CHECK_EQUAL( result.size(), 117000);
    std::ifstream myfile("tibra/tests/cpp_tests/results/element_classifier_cylinder.txt");
    std::string line;
    for( IndexType i = 0; i < result.size(); ++i){
        getline (myfile, line);
        BOOST_CHECK_EQUAL( result[i], std::stoi(line) );
        // Now test against flood fill solution.
        BOOST_CHECK_EQUAL( result[i], (*p_classification)[i]);
    }
    myfile.close();
}

BOOST_AUTO_TEST_CASE(CubeElementClassifierTest) {
    TIBRA_INFO << "Testing :: Test Classify Elements :: Cube with cavity" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

    Parameters params( {Component("lower_bound_xyz", PointType(-1.5, -1.5, -1.5)),
                        Component("upper_bound_xyz", PointType(1.5, 1.5, 1.5)),
                        Component("lower_bound_uvw", PointType(-1.5, -1.5, -1.5)),
                        Component("upper_bound_uvw", PointType(1.5, 1.5, 1.5)),
                        Component("number_of_elements", Vector3i(20, 20, 20)) });

    // Instatiate brep_operator
    BRepOperator brep_operator(triangle_mesh, params);
    Mapper mapper(params);

    std::vector<IndexType> result{};
    result.reserve(9261);
    for( IndexType i = 0; i < 8000; ++i ){
        const auto bounding_box = mapper.GetBoundingBoxXYZFromIndex(i);
        result.push_back( brep_operator.GetIntersectionState(bounding_box.first, bounding_box.second) );
    }

    // Get flood fill solution
    const auto p_classification = brep_operator.pGetElementClassifications();

    BOOST_CHECK_EQUAL( result.size(), 8000);
    std::ifstream myfile("tibra/tests/cpp_tests/results/element_classifier_cube.txt");
    std::string line;
    for( IndexType i = 0; i < result.size(); ++i){
        getline (myfile, line);
        BOOST_CHECK_EQUAL( result[i], std::stoi(line) );
        // Now test against flood fill solution.
        BOOST_CHECK_EQUAL( result[i], (*p_classification)[i]);
    }
    myfile.close();
}

BOOST_AUTO_TEST_CASE(ElephantElementClassifierTest) {
    TIBRA_INFO << "Testing :: Test Classify Elements :: Elephant" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/elephant.stl");

    Parameters params( {Component("lower_bound_xyz", PointType(-0.4, -0.6, -0.35)),
                        Component("upper_bound_xyz", PointType(0.4, 0.6, 0.35)),
                        Component("lower_bound_uvw", PointType(-0.4, -0.6, -0.35)),
                        Component("upper_bound_uvw", PointType(0.4, 0.6, 0.35)),
                        Component("number_of_elements", Vector3i(16, 24, 14)) });

    // Instatiate brep_operator
    BRepOperator brep_operator(triangle_mesh, params);
    Mapper mapper(params);

    std::vector<IndexType> result{};
    result.reserve(5376);
    for( IndexType i = 0; i < 5376; ++i ){
        const auto bounding_box = mapper.GetBoundingBoxXYZFromIndex(i);
        result.push_back( brep_operator.GetIntersectionState(bounding_box.first, bounding_box.second) );
    }

    // Get flood fill solution
    const auto p_classification = brep_operator.pGetElementClassifications();

    BOOST_CHECK_EQUAL( result.size(), 5376);
    std::ifstream myfile("tibra/tests/cpp_tests/results/element_classifier_elephant.txt");
    std::string line;
    for( IndexType i = 0; i < result.size(); ++i){
        getline (myfile, line);
        BOOST_CHECK_EQUAL( result[i], std::stoi(line) );
        // Now test against flood fill solution.
        BOOST_CHECK_EQUAL( result[i], (*p_classification)[i]);
    }
    myfile.close();
}

BOOST_AUTO_TEST_CASE(BunnyElementClassifierTest) {
    TIBRA_INFO << "Testing :: Test Classify Elements :: Bunny" << std::endl;

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/stanford_bunny.stl");

    Parameters params( {Component("lower_bound_xyz", PointType(-24, -43, 5)),
                        Component("upper_bound_xyz", PointType(85, 46, 115)),
                        Component("lower_bound_uvw", PointType(-24, -43, 5)),
                        Component("upper_bound_uvw", PointType(85, 46, 115)),
                        Component("number_of_elements", Vector3i(36, 30, 40)) });

    // Instatiate brep_operator
    BRepOperator brep_operator(triangle_mesh, params);
    Mapper mapper(params);

    std::vector<IndexType> result{};
    result.reserve(43200);
    for( IndexType i = 0; i < 43200; ++i){
        const auto bounding_box = mapper.GetBoundingBoxXYZFromIndex(i);
        result.push_back( brep_operator.GetIntersectionState(bounding_box.first, bounding_box.second ) );
    }

    // Get flood fill solution
    const auto p_classification = brep_operator.pGetElementClassifications();

    BOOST_CHECK_EQUAL( result.size(), 43200);
    std::ifstream myfile("tibra/tests/cpp_tests/results/element_classifier_bunny.txt");
    std::string line;
    for( IndexType i = 0; i < result.size(); ++i){
        getline (myfile, line);
        BOOST_CHECK_EQUAL( result[i], std::stoi(line) );
        // Now test against flood fill solution.
        BOOST_CHECK_EQUAL( result[i], (*p_classification)[i]);
    }
    myfile.close();
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra