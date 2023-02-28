// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// STL includes

//// Project includes
#include "containers/triangle_mesh.hpp"
#include "embedding/brep_operator.h"
#include "embedding/octree.h"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( Octree_TestSuite )

BOOST_AUTO_TEST_CASE(OctreeCubeTest1) {
    TIBRA_INFO << "Testing :: Test Octree :: Test Cube 1" << std::endl;
    /// Uniform refinement, each node has 8 children.

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

    Vector3d lower_bound = {0.0, 0.0, 0.0};
    Vector3d upper_bound = {2.0, 2.0, 2.0};
    Parameters params( {Component("lower_bound", lower_bound),
                        Component("upper_bound", upper_bound) });

    // Get trimmed domain.
    BRepOperator brep_operator(triangle_mesh, params);
    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain({-2.0, -2, -2},{-1.3, -1.3, -1.3});

    // Construct octree.
    Octree octree(p_trimmed_domain.get(), {-1.5, -1.5, -1.5},{-1.3, -1.3, -1.3} , params);

    // Refine octree to level 5.
    octree.Refine(4, 4);
    BOOST_CHECK_EQUAL(octree.NumberOfNodes(), 4681UL); // 8^0+8^1+8^2+8^3..+8^4
    BOOST_CHECK_EQUAL(octree.NumberOfLeafs(), 4096UL); // 8^4

    // Refine octree to level 6.
    octree.Refine(5, 5);
    BOOST_CHECK_EQUAL(octree.NumberOfNodes(), 37449UL); // 8^0+8^1+8^2+8^3..+8^5
    BOOST_CHECK_EQUAL(octree.NumberOfLeafs(), 32768UL); // 8^5

    Vector3i r_order(2, 3, 1);
    auto p_points = octree.pGetIntegrationPoints(r_order);
    BOOST_CHECK_EQUAL( p_points->size(), 786432 );
    double volume = 0.0;
    for( auto point : (*p_points)){
        volume += point.GetWeight()*8.0;
    }
    BOOST_CHECK_LT( std::abs(volume-0.008)/0.008, 1e-10);
} // End TouchingCubeTest1

BOOST_AUTO_TEST_CASE(OctreeCubeTest2) {
    TIBRA_INFO << "Testing :: Test Octree :: Test Cube 2" << std::endl;
    /// Refinement in only on trimmed nodes in one direction.

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

    Vector3d lower_bound = {-1.0, -1.0, -1.0};
    Vector3d upper_bound = {2.0, 2.0, 2.0};
    Parameters params( {Component("lower_bound", lower_bound),
                        Component("upper_bound", upper_bound) });

    // Get trimmed domain.
    BRepOperator brep_operator(triangle_mesh, params);
    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain({-2.0, -2, -2},{-1.3, -1.3, -1.3});

    // Construct octree.
    Octree octree(p_trimmed_domain.get(), {-1.50001, -1.49999, -1.49999},{-1.3, -1.3, -1.3} , params);

    // Refine octree to level 5.
    octree.Refine(0, 4);
    BOOST_CHECK_EQUAL(octree.NumberOfNodes(), 681UL);   // 1+4+4*4+4*4*4 + 4^1+4^2+4^3+4^4*2
    BOOST_CHECK_EQUAL(octree.NumberOfLeafs(), 596UL);   // 4^1+4^2+4^3+4^4*2

    Vector3i r_order(0, 0, 0);
    auto p_points = octree.pGetIntegrationPoints(r_order);

    BOOST_CHECK_EQUAL( p_points->size(), 596 );
    double volume = 0.0;
    for( auto point : (*p_points)){
        volume += point.GetWeight()*27.0;
    }
    BOOST_CHECK_LT( std::abs(volume-0.008)/0.008, 1e-4);
} // End OctreeCubeTest2

BOOST_AUTO_TEST_CASE(OctreeElephantTest) {
    TIBRA_INFO << "Testing :: Test Octree :: Test Elephant" << std::endl;
    // Compute volume of elephant through octree.

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/elephant.stl");

    Vector3d lower_bound = {0.0, 0.0, 0.0};
    Vector3d upper_bound = {5.0, 5.0, 5.0};
    Parameters params( {Component("lower_bound", lower_bound),
                        Component("upper_bound", upper_bound) });

    const double ref_volume = MeshUtilities::VolumeOMP(triangle_mesh);
    // Get trimmed domain.
    BRepOperator brep_operator(triangle_mesh, params);
    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain({-0.4, -0.6, -0.35},{0.4, 0.6, 0.35});

    // Construct octree.
    Octree octree(p_trimmed_domain.get(), {-0.4, -0.6, -0.35},{0.4, 0.6, 0.35} , params);

    // Refine octree to level 5.
    octree.Refine(0, 5);

    // Check if integration points contain same volume as ref_volume
    Vector3i r_order(2, 2, 2);
    auto p_points = octree.pGetIntegrationPoints(r_order);

    BOOST_CHECK_EQUAL( octree.NumberOfNodes(), 3887 );
    BOOST_CHECK_EQUAL( p_points->size(), 45186 );
    double volume = 0.0;
    for( auto point : (*p_points)){
        volume += point.GetWeight()*125.0;
    }
    BOOST_CHECK_LT( std::abs(volume - ref_volume) / ref_volume, 2e-4);

    // Refine inner levels -> volume must the same as before.
    octree.Refine(5, 5);
    auto p_points_2 = octree.pGetIntegrationPoints(r_order);
    BOOST_CHECK_EQUAL( p_points_2->size(), 60873);
    double volume_2 = 0.0;
    for( auto point : (*p_points_2)){
        volume_2 += point.GetWeight()*125.0;
    }
    BOOST_CHECK_LT( std::abs(volume - volume_2) / volume, 1e-10);
} // End OctreeBunny

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra

