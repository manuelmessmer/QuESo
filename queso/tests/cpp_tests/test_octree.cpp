//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//			    See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// STL includes

//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/brep_operator.h"
#include "queso/embedding/octree.h"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( Octree_TestSuite )

BOOST_AUTO_TEST_CASE(OctreeCubeTest1) {
    QuESo_INFO << "Testing :: Test Octree :: Test Cube 1" << std::endl;
    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    /// Uniform refinement, each node has 8 children.

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cube_with_cavity.stl");

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    // Get trimmed domain.
    BRepOperator brep_operator(triangle_mesh);
    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain({-2.0, -2, -2},{-1.3, -1.3, -1.3}, min_vol_ratio, min_num_triangles);

    // Construct octree.
    Octree<TrimmedDomain> octree(p_trimmed_domain.get(), MakeBox({-1.5, -1.5, -1.5},{-1.3, -1.3, -1.3}),
                                          MakeBox({-1.0, -1.0, -1.0},{1.0, 1.0, 1.0}) );

    // Refine octree to level 5.
    octree.Refine(4, 4);
    QuESo_CHECK_EQUAL(octree.NumberOfNodes(), 4681UL); // 8^0+8^1+8^2+8^3..+8^4
    QuESo_CHECK_EQUAL(octree.NumberOfLeafs(), 4096UL); // 8^4

    // Refine octree to level 6.
    octree.Refine(5, 5);
    QuESo_CHECK_EQUAL(octree.NumberOfNodes(), 37449UL); // 8^0+8^1+8^2+8^3..+8^5
    QuESo_CHECK_EQUAL(octree.NumberOfLeafs(), 32768UL); // 8^5

    Vector3i r_order{2, 3, 1};
    auto p_points = octree.pGetIntegrationPoints<ElementType>(r_order);
    QuESo_CHECK_EQUAL( p_points->size(), 786432 );
    double volume = 0.0;
    for( auto point : (*p_points)){
        volume += point.Weight();
    }
    QuESo_CHECK_LT( std::abs(volume-8.0)/8.0, 1e-10);
} // End TouchingCubeTest1

BOOST_AUTO_TEST_CASE(OctreeCubeTest2) {
    QuESo_INFO << "Testing :: Test Octree :: Test Cube 2" << std::endl;

    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;
    /// Refinement in only on trimmed nodes in one direction.

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cube_with_cavity.stl");

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    // Get trimmed domain.
    BRepOperator brep_operator(triangle_mesh);
    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain({-2.0, -2, -2},{-1.3, -1.3, -1.3}, min_vol_ratio, min_num_triangles);

    // Construct octree.
    Octree<TrimmedDomain> octree(p_trimmed_domain.get(), MakeBox({-1.50001, -1.49999, -1.49999},{-1.3, -1.3, -1.3}),
                                          MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));

    // Refine octree to level 5.
    octree.Refine(0, 4);
    QuESo_CHECK_EQUAL(octree.NumberOfNodes(), 681UL);   // 1+4+4*4+4*4*4 + 4^1+4^2+4^3+4^4*2
    QuESo_CHECK_EQUAL(octree.NumberOfLeafs(), 596UL);   // 4^1+4^2+4^3+4^4*2

    Vector3i r_order{0, 0, 0};
    auto p_points = octree.pGetIntegrationPoints<ElementType>(r_order);

    QuESo_CHECK_EQUAL( p_points->size(), 596 );
    double volume = 0.0;
    for( auto point : (*p_points)){
        volume += point.Weight();
    }
    QuESo_CHECK_LT( std::abs(volume-1.0)/1.0, 1e-4);
} // End OctreeCubeTest2

BOOST_AUTO_TEST_CASE(OctreeElephantTest) {
    QuESo_INFO << "Testing :: Test Octree :: Test Elephant" << std::endl;
    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    // Compute volume of elephant through octree.

    // Read mesh from STL file
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/elephant.stl");

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 500;

    const double ref_volume = MeshUtilities::VolumeOMP(triangle_mesh);
    // Get trimmed domain.
    BRepOperator brep_operator(triangle_mesh);
    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain({-0.4, -0.6, -0.35},{0.4, 0.6, 0.35}, min_vol_ratio, min_num_triangles);

    // Construct octree.
    Octree<TrimmedDomain> octree(p_trimmed_domain.get(), MakeBox({-0.4, -0.6, -0.35},{0.4, 0.6, 0.35}),
                                          MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}));

    // Refine octree to level 5.
    octree.Refine(0, 5);

    // Check if integration points contain same volume as ref_volume
    Vector3i r_order{2, 2, 2};
    auto p_points = octree.pGetIntegrationPoints<ElementType>(r_order);

    QuESo_CHECK_EQUAL( octree.NumberOfNodes(), 3887 );
    QuESo_CHECK_EQUAL( p_points->size(), 45186 );
    double volume = 0.0;
    for( auto point : (*p_points)){
        volume += point.Weight()*(0.8*1.2*0.7) / 8.0;;
    }
    QuESo_CHECK_LT( std::abs(volume - ref_volume) / ref_volume, 2e-4);

    // Refine inner levels -> volume must the same as before.
    octree.Refine(5, 5);
    auto p_points_2 = octree.pGetIntegrationPoints<ElementType>(r_order);
    QuESo_CHECK_EQUAL( p_points_2->size(), 60873);
    double volume_2 = 0.0;
    for( auto point : (*p_points_2)){
        volume_2 += point.Weight()*(0.8*1.2*0.7) / 8.0;
    }
    QuESo_CHECK_LT( std::abs(volume - volume_2) / volume, 1e-10);
} // End OctreeBunny

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

