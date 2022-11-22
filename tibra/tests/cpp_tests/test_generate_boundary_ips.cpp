// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <chrono>

#include "geometries/triangle_mesh.h"
#include "io/io_utilities.h"
#include "embedding/brep_operator.h"

#include "cgal_wrapper/cgal_brep_operator.h"
#include "quadrature/moment_fitting_utilities.h"



#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/STL.h>
#include <boost/numeric/ublas/matrix.hpp>

namespace Testing{

BOOST_AUTO_TEST_SUITE( GenerateBoundaryIPsTestSuite )

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsTest1) {
    std::cout << "Testing :: Prototype :: Generate Boundary Integration Points :: Elephant" << std::endl;
    /// Test of prototype Functions.
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3;
    typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
    typedef boost::numeric::ublas::vector<double> VectorType;

    SurfaceMeshType mPolyhedron;
    CGAL::IO::read_STL("tibra/tests/cpp_tests/data/elephant.stl", mPolyhedron);

    TriangleMesh::Vector3d lower_bound = {0.1, 0.1, -0.1};
    TriangleMesh::Vector3d upper_bound = {0.2, 0.2, 0.1};

    std::array<int, 3> number_of_elements = {1, 1, 1};
    std::array<int, 3> order = {2, 2, 2};

    int point_distribution_factor = 3;
    double initial_triangle_edge_length = 1;
    int minimum_number_of_triangles = 10000;
    double moment_fitting_residual = 1e-8;
    std::string integration_method = "Gauss";
    int echo_level = 0;

    Parameters param(lower_bound, upper_bound, number_of_elements, order, initial_triangle_edge_length,
        minimum_number_of_triangles, moment_fitting_residual, point_distribution_factor, integration_method, echo_level);

    Element element(1, lower_bound, upper_bound, param);

    // auto status = cgal::BRepOperator::ComputeIntersectionMesh( mPolyhedron, *p_cube, element, param);

    // VectorType constant_terms_cgal{};
    // cgal::ConstantTerms::Compute(element, constant_terms_cgal, param);


    // //Read mesh from STL file
    // TriangleMesh triangle_mesh{};
    // IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/elephant.stl");

    // // Build brep_operator
    // BRepOperator brep_operator(triangle_mesh);
    // auto p_points = brep_operator.GetBoundaryIps(lower_bound, upper_bound);

    // VectorType constant_terms{};
    // ConstantTerms::Compute(p_points, element, constant_terms, param);


    // // Check
    // for( int i = 0; i < constant_terms_cgal.size(); ++i){
    //     double ralative_error = (constant_terms[i] - constant_terms_cgal[i])/constant_terms_cgal[i];
    //     //std::cout << ralative_error << "\t " << constant_terms[i] << ", " << constant_terms_cgal[i] << std::endl;
    //     BOOST_CHECK_SMALL(ralative_error, 0.03);

    // }

    // auto cube = Modeler::make_cube_3(lower_bound, upper_bound);
    // // Clip mesh.
    // auto test = brep_operator.GetBoundaryIps(lower_bound, upper_bound);
    // auto p_clipped_mesh = brep_operator.ClipTriangleMesh(lower_bound, upper_bound);

    //IO::WriteMeshToVTK(*cube, "cube.vtk", true);
    //IO::WriteMeshToSTL(*p_clipped_mesh, "test.stl", true);

}

BOOST_AUTO_TEST_SUITE_END()

} // End GenerateBoundaryIPsTestSuite
