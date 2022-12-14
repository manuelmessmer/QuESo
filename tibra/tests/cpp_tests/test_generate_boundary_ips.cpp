// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <chrono>

#include "containers/triangle_mesh.h"
#include "io/io_utilities.h"
#include "embedding/brep_operator.h"
#include "quadrature/moment_fitting_utilities.h"

#include "cgal_wrapper/cgal_brep_operator.h"
#include "cgal_wrapper/cgal_io_utilities.h"
#include "cgal_wrapper/cgal_cuboid_modeler.h"
// #include "quadrature/moment_fitting_utilities.h"


namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( GenerateBoundaryIPsTestElephantSuite )

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsTest1) {
    std::cout << "Testing :: Prototype :: Generate Boundary Integration Points :: Elephant" << std::endl;

    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};
    Vector3d lower_bound = {-0.4, -0.6, -0.35 };
    Vector3d upper_bound = {0.4, 0.6, 0.35 };
    Parameters parameters( {Component("lower_bound", lower_bound),
                            Component("upper_bound", upper_bound),
                            Component("min_num_boundary_triangles", 200UL),
                            Component("number_of_elements", number_of_elements),
                            Component("polynomial_order", Vector3i(2,2,2) ) } );

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/elephant.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh,parameters);

    const double delta_x = 0.1;
    const double delta_y = 0.1;
    const double delta_z = 0.1;

    std::ifstream file("tibra/tests/cpp_tests/results/surface_integral_elephant.txt");
    std::string line{};

    auto start_time = std::chrono::high_resolution_clock::now();
    IndexType number_trimmed_elements = 0;

    for(double x = -0.4; x <= 0.4; x += delta_x){
        for(double y = -0.6; y <= 0.6; y += delta_y){
            for(double z = -0.35; z <= 0.35; z += delta_x){
                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto local_lower_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_lower_bound, lower_bound, upper_bound);
                auto local_upper_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_upper_bound, lower_bound, upper_bound);

                Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);
                if( brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound) == BRepOperator::Trimmed){
                    // Get trimmed domain
                    auto p_trimmed_domain = brep_operator.GetTrimmedDomain(local_lower_bound, local_upper_bound);
                    // Get boundary integration points
                    auto p_boundary_ips = p_trimmed_domain->pGetBoundaryIps();

                    VectorType constant_terms{};
                    MomentFitting::ComputeConstantTerms(element, p_boundary_ips, constant_terms, parameters);

                    // // Read and ignore header
                    getline(file, line);

                    double surface_area = 0.0;
                    for( auto& ip : *p_boundary_ips){
                        surface_area += ip.GetWeight();
                    }

                    // Read ref surface area
                    getline(file, line);
                    double ref_surface_area = std::stod(line);
                    BOOST_CHECK_LT(  std::abs(surface_area-ref_surface_area)/std::abs(ref_surface_area), 1e-10);

                    double error = 0.0;
                    double norm = 0.0;
                    for( auto& value : constant_terms ){
                        getline(file, line);
                        double ref_value = std::stod(line);
                        error += std::abs(value-ref_value);
                        norm += std::abs(ref_value);
                    }

                    BOOST_CHECK_LT( error/norm, 1e-6);
                    number_trimmed_elements++;
                }
            }
        }
    }
    file.close();
    BOOST_CHECK_EQUAL(number_trimmed_elements, 166);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;
    std::cout << "TIBRA :: Elapsed Time: " << elapsed_time.count() << std::endl;
}

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsBunnyTest) {
    std::cout << "Testing :: Prototype :: Generate Boundary Integration Points :: Bunny" << std::endl;

    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};
    Vector3d lower_bound = {-0.4, -0.6, -0.35 };
    Vector3d upper_bound = {0.4, 0.6, 0.35 };
    Parameters parameters( {Component("lower_bound", lower_bound),
                            Component("upper_bound", upper_bound),
                            Component("number_of_elements", number_of_elements),
                            Component("min_num_boundary_triangles", 100UL),
                            Component("polynomial_order", Vector3i(2,2,2) ) } );
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/stanford_bunny.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh, parameters);

    const double delta_x = 15;
    const double delta_y = 15;
    const double delta_z = 15;

    //auto start_time = std::chrono::high_resolution_clock::now();

    std::ifstream file("tibra/tests/cpp_tests/results/surface_integral_bunny.txt");
    std::string line{};

    IndexType number_trimmed_elements = 0;
    for(double x = -24; x <= 85; x += delta_x){
        for(double y = -43; y <= 46; y += delta_y){
            for(double z = 5; z <= 115; z += delta_z){
                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto local_lower_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_lower_bound, lower_bound, upper_bound);
                auto local_upper_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_upper_bound, lower_bound, upper_bound);
                Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);
                if( brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound) == BRepOperator::Trimmed){
                    // Get Trimmed domain
                    auto p_trimmed_domain = brep_operator.GetTrimmedDomain(local_lower_bound, local_upper_bound);
                    // Get boundary integration points
                    auto p_boundary_ips = p_trimmed_domain->pGetBoundaryIps();
                    // Read and ignore header
                    getline(file, line);

                    double area = 0.0;
                    for( const auto& ip : *p_boundary_ips){
                        area += ip.GetWeight();
                    }
                    getline(file, line);
                    const double ref_area = stod(line);
                    BOOST_CHECK_LT( std::abs(area-ref_area)/std::abs(ref_area), 1e-10);

                    VectorType constant_terms{};
                    MomentFitting::ComputeConstantTerms(element, p_boundary_ips, constant_terms, parameters);

                    double error = 0.0;
                    double norm = 0.0;
                    for( const auto& value : constant_terms ){
                        getline(file, line);
                        double ref_value = std::stod(line);
                        error += std::abs(value-ref_value);
                        norm += std::abs(ref_value);
                    }

                    BOOST_CHECK_LT( error/norm, 1e-6);
                    number_trimmed_elements++;
                }
            }
        }
    }
    file.close();
    // auto end_time = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed_time = end_time - start_time;
    // std::cout << "TIBRA :: Elapsed Time: " << elapsed_time.count() << std::endl;
    // std::cout << "number_trimmed_elements: " << number_trimmed_elements << std::endl;
}



// BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsTest2) {
//     std::cout << "Testing :: Prototype :: Generate Boundary Integration Points :: Bunny" << std::endl;

//     typedef boost::numeric::ublas::vector<double> VectorType;

//     Vector3i number_of_elements = {1, 1, 1};
//     Vector3d lower_bound = {-0.4, -0.6, -0.35 };
//     Vector3d upper_bound = {0.4, 0.6, 0.35 };
//     Parameters parameters( {Component("lower_bound", lower_bound),
//                             Component("upper_bound", upper_bound),
//                             Component("number_of_elements", number_of_elements),
//                             Component("min_num_boundary_triangles", 100UL),
//                             Component("polynomial_order", Vector3i(2,2,2) ) } );
//     TriangleMesh triangle_mesh{};
//     IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/stanford_bunny.stl");

//     // Build brep_operator
//     cgal::CGALBRepOperator cgal_operator_cgal(triangle_mesh, parameters);
//     BRepOperator brep_operator(triangle_mesh, parameters);

//     const double delta_x = 15;
//     const double delta_y = 15;
//     const double delta_z = 15;

//     auto start_time = std::chrono::high_resolution_clock::now();

//     std::ofstream file;
//     file.open("results.txt");

//     IndexType number_trimmed_elements = 0;
//     for(double x = -24; x <= 85; x += delta_x){
//         for(double y = -43; y <= 46; y += delta_y){
//             for(double z = 5; z <= 115; z += delta_z){
//                 Vector3d local_lower_bound = {x, y, z};
//                 Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

//                 auto local_lower_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_lower_bound, lower_bound, upper_bound);
//                 auto local_upper_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_upper_bound, lower_bound, upper_bound);
//                 Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);
//                 if( brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound) == BRepOperator::Trimmed){
//                     auto p_trimmed_domain_cgal = cgal_operator_cgal.GetTrimmedDomain(local_lower_bound, local_upper_bound, parameters);
//                     if( p_trimmed_domain_cgal ){
//                         auto p_boundary_ips_cgal = p_trimmed_domain_cgal->pGetBoundaryIps();

//                         VectorType constant_terms_cgal{}; // if( std::abs(distance) < 1e-10 && !Positive){
//                         MomentFitting::ComputeConstantTerms(element, p_boundary_ips_cgal, constant_terms_cgal, parameters);

//                         auto p_trimmed_domain = brep_operator.GetTrimmedDomain(local_lower_bound, local_upper_bound);
//                         const auto& r_clipped_mesh = p_trimmed_domain->GetTriangleMesh();
//                         //IO::WriteMeshToSTL(r_clipped_mesh, "clipped_mesh_before.stl", true);
//                         auto p_boundary_ips = p_trimmed_domain->pGetBoundaryIps();
//                         IO::WriteMeshToSTL(r_clipped_mesh, "clipped_mesh.stl", true);

//                         const auto& r_clipped_mesh_cgal = p_trimmed_domain_cgal->GetTriangleMesh();
//                         IO::WriteMeshToSTL(r_clipped_mesh_cgal, "clipped_mesh_cgal.stl", true);

//                         IO::WritePointsToVTK(*p_boundary_ips, "boundary_ips.vtk", true);
//                         IO::WritePointsToVTK(*p_boundary_ips_cgal, "boundary_ips_cgal.vtk", true);

//                         file << "E: " << number_trimmed_elements << '\n';
//                         file << std::setprecision(16);
//                         double area = 0.0;
//                         for( auto& ip : *p_boundary_ips){
//                             area += ip.GetWeight();
//                         }
//                         double area_cgal = 0.0;
//                         for( auto& ip : *p_boundary_ips_cgal){
//                             area_cgal += ip.GetWeight();
//                         }
//                         double area_error = std::abs(area-area_cgal)/std::abs(area_cgal);
//                         file << area << '\n';
//                         //std::cout << "area_error: " << area_error << std::endl;
//                         BOOST_CHECK_LT(area_error, 1e-8);

//                         VectorType constant_terms{};
//                         MomentFitting::ComputeConstantTerms(element, p_boundary_ips, constant_terms, parameters);

//                         //std::cout << "Volume" << ": " << constant_terms_cgal[0] << std::endl;
//                         double volume_error = std::abs(constant_terms[0] - constant_terms_cgal[0] )/ std::abs(constant_terms_cgal[0]);

//                         if( constant_terms[0] > 1e-3 ){
//                             std::cout << volume_error << std::endl;
//                             BOOST_CHECK_LT(volume_error, 2e-5);

//                         }
//                         double val_error = 0.0;
//                         double norm = 0.0;
//                         for( IndexType i = 0; i < constant_terms.size(); ++i){
//                             file << constant_terms[i] << std::endl;
//                             val_error += std::abs(constant_terms[i]-constant_terms_cgal[i]);
//                             norm += std::abs(constant_terms_cgal[i]);
//                         }
//                         std::cout << "Val error: " << val_error/norm << std::endl;
//                         BOOST_CHECK_LT(val_error/norm, 1e-4);
//                         // auto cube = cgal::CuboidModeler::MakeCuboid(local_lower_bound, local_upper_bound);
//                         // cgal::IO::WriteMeshToVTK(*cube, "cube.vtk", true);

//                         number_trimmed_elements++;

//                         //std::cout << "number_trimmed_elements: " << number_trimmed_elements << std::endl;
//                     }


//                 }
//             }
//         }
//     }
//     file.close();
//     auto end_time = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_time = end_time - start_time;
//     std::cout << "TIBRA :: Elapsed Time: " << elapsed_time.count() << std::endl;
//     std::cout << "number_trimmed_elements: " << number_trimmed_elements << std::endl;

// }

// BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsTest2) {
//     std::cout << "Testing :: Prototype :: Generate Boundary Integration Points :: Bunny" << std::endl;
//     /// Test of prototype Functions.
//     // typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//     // typedef K::Point_3 Point_3;
//     // typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
//     typedef boost::numeric::ublas::vector<double> VectorType;

//     // SurfaceMeshType mPolyhedron;
//     // CGAL::IO::read_STL("tibra/tests/cpp_tests/data/elephant.stl", mPolyhedron);

//     Vector3i number_of_elements = {1, 1, 1};

//     Vector3d lower_bound = {-0.4, -0.6, -0.35 };
//     Vector3d upper_bound = {0.4, 0.6, 0.35 };
//     Parameters parameters( {Component("lower_bound", lower_bound),
//                             Component("upper_bound", upper_bound),
//                             Component("number_of_elements", number_of_elements),
//                             Component("min_num_boundary_triangles", 100UL),
//                             Component("polynomial_order", Vector3i(2,2,2) ) } );
//     TriangleMesh triangle_mesh{};
//     IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/stanford_bunny.stl");

//     // Build brep_operator
//     cgal::CGALBRepOperator cgal_operator_cgal(triangle_mesh);
//     BRepOperator brep_operator(triangle_mesh);

//     const double delta_x = 3;
//     const double delta_y = 3;
//     const double delta_z = 3;
//     auto start_time = std::chrono::high_resolution_clock::now();
//     IndexType number_trimmed_elements = 0;
//     for(double x = -24; x <= 85; x += delta_x){
//         for(double y = -43; y <= 46; y += delta_y){
//             for(double z = 5; z <= 115; z += delta_z){
//                 Vector3d local_lower_bound = {x, y, z};
//                 Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

//                 auto local_lower_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_lower_bound, lower_bound, upper_bound);
//                 auto local_upper_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_upper_bound, lower_bound, upper_bound);
//                 Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);
//                 if( brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound) == BRepOperator::Trimmed){



//                         auto p_trimmed_domain = brep_operator.GetTrimmedDomain(local_lower_bound, local_upper_bound);
//                         auto p_boundary_ips = p_trimmed_domain->pGetBoundaryIps();


//                         VectorType constant_terms{};
//                         //MomentFitting::ComputeConstantTerms(element, p_boundary_ips, constant_terms, parameters);
//                         number_trimmed_elements++;
//                         //std::cout << "number_trimmed_elements: " << number_trimmed_elements << std::endl;
//                     //}
//                     // if( number_trimmed_elements == 6 ){
//                     //     throw std::runtime_error("waf");
//                     //}

//                 }
//             }
//         }
//     }
//     std::cout << "number_trimmed_elements: " << number_trimmed_elements << std::endl;
//     auto end_time = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_time = end_time - start_time;
//     std::cout << "TIBRA :: Elapsed Time: " << elapsed_time.number_trimmed_elements() << std::endl;

// }

// BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsTest3) {
//     std::cout << "Testing :: Prototype :: Generate Boundary Integration Points :: Cube" << std::endl;
//     /// Test of prototype Functions.
//     // typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//     // typedef K::Point_3 Point_3;
//     // typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
//     typedef boost::numeric::ublas::vector<double> VectorType;

//     // SurfaceMeshType mPolyhedron;
//     // CGAL::IO::read_STL("tibra/tests/cpp_tests/data/elephant.stl", mPolyhedron);

//     Vector3i number_of_elements = {1, 1, 1};
//     // Vector3i order = {2, 2, 2};

//     // int point_distribution_factor = 3;
//     // double initial_triangle_edge_length = 1;
//     // int minimum_number_of_triangles = 10000;
//     // double moment_fitting_residual = 1e-8;
//     // std::string integration_method = "Gauss";
//     // int echo_level = 0;
//     Vector3d lower_bound = {-2, -2, -2 };
//     Vector3d upper_bound = {-1, -1, -1 };
//     Parameters parameters( {Component("lower_bound", lower_bound),
//                             Component("upper_bound", upper_bound),
//                             Component("number_of_elements", number_of_elements),
//                             Component("min_num_boundary_triangles", 50000UL),
//                             Component("polynomial_order", Vector3i(2,2,2) ) } );
//     TriangleMesh triangle_mesh{};
//     IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

//     // Build brep_operator
//     cgal::CGALBRepOperator cgal_operator_cgal(triangle_mesh);
//     BRepOperator brep_operator(triangle_mesh);

//     const double delta_x = 0.2;
//     const double delta_y = 0.2;
//     const double delta_z = 0.2;


//     Vector3d local_lower_bound = {-2, -2, -2};
//     Vector3d local_upper_bound = {-1, -1, -1};

//     auto local_lower_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_lower_bound, lower_bound, upper_bound);
//     auto local_upper_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_upper_bound, lower_bound, upper_bound);
//     Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);



//     auto p_trimmed_domain_cgal = cgal_operator_cgal.GetTrimmedDomain(local_lower_bound, local_upper_bound, parameters);
//     if( !p_trimmed_domain_cgal ){
//         std::cout << "p_trimmed_domain_cgal  is nullptr. \n";
//     }
//     auto p_boundary_ips_cgal = p_trimmed_domain_cgal->pGetBoundaryIps();

//     VectorType constant_terms_cgal{};
//     MomentFitting::ComputeConstantTerms(element, p_boundary_ips_cgal, constant_terms_cgal, parameters);

//     auto p_trimmed_domain = brep_operator.GetTrimmedDomain(local_lower_bound, local_upper_bound);
//     auto p_boundary_ips = p_trimmed_domain->pGetBoundaryIps();
//     const auto& r_clipped_mesh = p_trimmed_domain->GetTriangleMesh();
//     IO::WriteMeshToSTL(r_clipped_mesh, "clipped_mesh.stl", true);

//     const auto& r_clipped_mesh_cgal = p_trimmed_domain_cgal->GetTriangleMesh();
//     IO::WriteMeshToSTL(r_clipped_mesh_cgal, "clipped_mesh_cgal.stl", true);

//     IO::WritePointsToVTK(*p_boundary_ips, "boundary_ips.vtk", true);
//     IO::WritePointsToVTK(*p_boundary_ips_cgal, "boundary_ips_cgal.vtk", true);

//     VectorType constant_terms{};
//     MomentFitting::ComputeConstantTerms(element, p_boundary_ips, constant_terms, parameters);

//     double error_norm = 0.0;
//     for( IndexType i = 0; i < constant_terms.size(); ++i){
//         error_norm += std::abs( constant_terms[i] - constant_terms_cgal[i] );
//         std::cout << i << ": " << std::abs( constant_terms[i] - constant_terms_cgal[i] ) << std::endl;
//         std::cout << i << ": " << constant_terms[i] << ", Ref: " << constant_terms_cgal[i] << std::endl;
//     }
//     error_norm /= constant_terms.size();
//     std::cout << "error_norm: " <<  error_norm << std::endl;

//     auto cube = cgal::CuboidModeler::MakeCuboid(local_lower_bound, local_upper_bound);
//     cgal::IO::WriteMeshToVTK(*cube, "cube.vtk", true);
//     BOOST_CHECK_LT(error_norm, 1e-4);
//     //number_trimmed_elements++;
//     // if( error_norm > 1e-8 ){
//     //     std::cout << "number_trimmed_elements: " << number_trimmed_elements << std::endl;
//     //     throw std::runtime_error("waf");
//     // }
//     //std::cout << "number_trimmed_elements: " << number_trimmed_elements << std::endl;
//     // if( number_trimmed_elements == 5 ){
//     //     throw std::runtime_error("waf");
//     // }
// }

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra
