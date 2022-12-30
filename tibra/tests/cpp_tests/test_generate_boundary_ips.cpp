// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <chrono>

#include "containers/triangle_mesh.h"
#include "io/io_utilities.h"
#include "embedding/brep_operator.h"
#include "quadrature/moment_fitting_utilities.h"


namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( GenerateBoundaryIPsTestSuite )

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsElephantTest) {
    std::cout << "Testing :: Generate Boundary Integration Points :: Elephant" << std::endl;

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

    //auto start_time = std::chrono::high_resolution_clock::now();
    IndexType number_trimmed_elements = 0;

    for(double x = lower_bound[0]; x <= upper_bound[0]; x += delta_x){
        for(double y = lower_bound[1]; y <= upper_bound[1]; y += delta_y){
            for(double z = lower_bound[2]; z <= upper_bound[2]; z += delta_z){
                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto local_lower_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_lower_bound, lower_bound, upper_bound);
                auto local_upper_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_upper_bound, lower_bound, upper_bound);

                Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);
                if( brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound) == BRepOperatorBase::Trimmed){
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

    // auto end_time = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed_time = end_time - start_time;
    // std::cout << "TIBRA :: Elapsed Time: " << elapsed_time.count() << std::endl;
}

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsBunnyTest) {
    std::cout << "Testing :: Generate Boundary Integration Points :: Bunny" << std::endl;

    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};
    Vector3d lower_bound = {-24.0, -43.0, 5.0 };
    Vector3d upper_bound = {-5.0, 46.0, 115 };
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
    for(double x = lower_bound[0]; x <= upper_bound[0]; x += delta_x){
        for(double y = lower_bound[1]; y <= upper_bound[1]; y += delta_y){
            for(double z = lower_bound[2]; z <= upper_bound[2]; z += delta_z){
                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto local_lower_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_lower_bound, lower_bound, upper_bound);
                auto local_upper_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_upper_bound, lower_bound, upper_bound);
                Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);
                if( brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound) == BRepOperatorBase::Trimmed){
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

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsCylinderTest) {
    std::cout << "Testing :: Generate Boundary Integration Points :: Cylinder" << std::endl;

    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};

    Vector3d lower_bound = {-1.5, -1.5, -1 };
    Vector3d upper_bound = {1.5, 1.5, 12 };
    Parameters parameters( {Component("lower_bound", lower_bound),
                            Component("upper_bound", upper_bound),
                            Component("number_of_elements", number_of_elements),
                            Component("min_num_boundary_triangles", 100UL),
                            Component("polynomial_order", Vector3i(2,2,2) ) } );

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh, parameters);

    const double delta_x = 1;
    const double delta_y = 1;
    const double delta_z = 1;

    std::ifstream file("tibra/tests/cpp_tests/results/surface_integral_cylinder.txt");
    std::string line{};

    IndexType number_trimmed_elements = 0;
    for(double x = lower_bound[0]; x <= upper_bound[0]; x += delta_x){
        for(double y = lower_bound[1]; y <= upper_bound[1]; y += delta_y){
            for(double z = lower_bound[2]; z <= upper_bound[2]; z += delta_z){
                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto local_lower_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_lower_bound, lower_bound, upper_bound);
                auto local_upper_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_upper_bound, lower_bound, upper_bound);
                Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);

                if( brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound) == BRepOperatorBase::Trimmed){
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
                    for( IndexType i = 0; i < constant_terms.size(); ++i ){
                        double value = constant_terms[i];
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

    BOOST_CHECK_EQUAL( number_trimmed_elements, 80);
}


void RunCube(const PointType rDelta, const PointType rLowerBound, const PointType rUpperBound,
    const PointType Perturbation, const std::string rFilename ){
    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};
    Parameters parameters( {Component("lower_bound", rLowerBound),
                            Component("upper_bound", rUpperBound),
                            Component("number_of_elements", number_of_elements),
                            Component("min_num_boundary_triangles", 100UL),
                            Component("polynomial_order", Vector3i(2,2,2) ) } );

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

    auto& vertices = triangle_mesh.GetVertices();
    for( auto& v : vertices ){
        v[0] += Perturbation[0];
        v[1] += Perturbation[1];
        v[2] += Perturbation[2];
    }

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh, parameters);

    const double delta_x = rDelta[0];
    const double delta_y = rDelta[1];
    const double delta_z = rDelta[2];

    std::string filename = "tibra/tests/cpp_tests/results/b_ips_cube/" + rFilename + ".txt";
    std::ifstream file(filename);
    std::string line{};

    IndexType number_trimmed_elements = 0;
    for(double x = rLowerBound[0]; x <= rUpperBound[0]; x += delta_x){
        for(double y = rLowerBound[1]; y <= rUpperBound[1]; y += delta_y){
            for(double z = rLowerBound[2]; z <= rUpperBound[2]; z += delta_z){
                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto local_lower_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_lower_bound, rLowerBound, rUpperBound);
                auto local_upper_bound_param = MappingUtilities::FromGlobalToLocalSpace(local_upper_bound, rLowerBound, rUpperBound);
                Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);

                // Get Trimmed domain
                auto p_trimmed_domain = brep_operator.GetTrimmedDomain(local_lower_bound, local_upper_bound);
                if( p_trimmed_domain ){

                    // Get boundary integration points
                    auto p_boundary_ips = p_trimmed_domain->pGetBoundaryIps();

                    // Read and ignore header
                    getline(file, line);
                    double area = 0.0;
                    for( const auto& ip : *p_boundary_ips){
                        area += ip.GetWeight();
                    }
                    getline(file, line);
                    double ref_area = std::stod(line);

                    double area_error = std::abs(area-ref_area)/std::abs(ref_area);
                    BOOST_CHECK_LT( area_error, 5e-8 );

                    VectorType constant_terms{};
                    MomentFitting::ComputeConstantTerms(element, p_boundary_ips, constant_terms, parameters);

                    double error = 0.0;
                    double norm = 0.0;

                    double value_0 = 1.0;
                    for( IndexType i = 0; i < constant_terms.size(); ++i){
                        getline(file,line);
                        double value = std::stod(line);
                        if( i == 0 ){
                            value_0 = value;
                        }
                        error += std::abs( constant_terms[i] - value );
                        norm += std::abs( value );
                    }
                    if( std::abs(value_0) > 1e-5 ){
                        BOOST_CHECK_LT(error/norm, 1e-6);
                    }

                    number_trimmed_elements++;
                }
            }
        }
    }

    getline(file, line);
    IndexType value = static_cast<IndexType>(std::stoi(line));
    BOOST_CHECK_EQUAL(value, number_trimmed_elements);
    file.close();
}

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsCube1Test) {
    std::cout << "Testing :: Generate Boundary Integration Points :: Cube 1" << std::endl;

    std::vector<double> perturbations = { 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 };
    std::vector<std::string> filenames = { "positive_x_e6",
                                           "positive_x_e7",
                                           "positive_x_e8",
                                           "positive_x_e9",
                                           "positive_x_e10",
                                           "positive_x_e11",
                                           "positive_x_e12" };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-1.5, -1.5, -1.5};
        PointType upper_bound = {1.5, 1.5, 1.5};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {perturbations[i], 0.0 , 0.0};
        RunCube(delta, lower_bound, upper_bound, perturbation, filenames[i]);
    }
}

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsCube2Test) {
    std::cout << "Testing :: Generate Boundary Integration Points :: Cube 2" << std::endl;

    std::vector<double> perturbations = { 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 };
    std::vector<std::string> filenames = { "positive_y_e6",
                                           "positive_y_e7",
                                           "positive_y_e8",
                                           "positive_y_e9",
                                           "positive_y_e10",
                                           "positive_y_e11",
                                           "positive_y_e12" };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-1.5, -1.5, -1.5};
        PointType upper_bound = {1.5, 1.5, 1.5};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {0.0, perturbations[i], 0.0};
        RunCube(delta, lower_bound, upper_bound, perturbation, filenames[i]);
    }
}

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsCube3Test) {
    std::cout << "Testing :: Generate Boundary Integration Points :: Cube 3" << std::endl;

    std::vector<double> perturbations = { 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 };
    std::vector<std::string> filenames = { "positive_z_e6",
                                           "positive_z_e7",
                                           "positive_z_e8",
                                           "positive_z_e9",
                                           "positive_z_e10",
                                           "positive_z_e11",
                                           "positive_z_e12" };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-1.5, -1.5, -1.5};
        PointType upper_bound = {1.5, 1.5, 1.5};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {0.0, 0.0, perturbations[i]};
        RunCube(delta, lower_bound, upper_bound, perturbation, filenames[i]);
    }
}

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsCube4Test) {
    std::cout << "Testing :: Generate Boundary Integration Points :: Cube 4" << std::endl;

    std::vector<double> perturbations = { 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 };
    std::vector<std::string> filenames = { "negative_x_e6",
                                           "negative_x_e7",
                                           "negative_x_e8",
                                           "negative_x_e9",
                                           "negative_x_e10",
                                           "negative_x_e11",
                                           "negative_x_e12" };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-1.5, -1.5, -1.5};
        PointType upper_bound = {1.5, 1.5, 1.5};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {-perturbations[i], 0.0 , 0.0};
        RunCube(delta, lower_bound, upper_bound, perturbation, filenames[i]);
    }
}

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsCube5Test) {
    std::cout << "Testing :: Generate Boundary Integration Points :: Cube 5" << std::endl;

    std::vector<double> perturbations = { 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 };
    std::vector<std::string> filenames = { "negative_y_e6",
                                           "negative_y_e7",
                                           "negative_y_e8",
                                           "negative_y_e9",
                                           "negative_y_e10",
                                           "negative_y_e11",
                                           "negative_y_e12" };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-1.5, -1.5, -1.5};
        PointType upper_bound = {1.5, 1.5, 1.5};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {0.0, -perturbations[i], 0.0};
        RunCube(delta, lower_bound, upper_bound, perturbation, filenames[i]);
    }
}

BOOST_AUTO_TEST_CASE(GenerateBoundaryIPsCube6Test) {
    std::cout << "Testing :: Generate Boundary Integration Points :: Cube 6" << std::endl;

    std::vector<double> perturbations = { 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 };
    std::vector<std::string> filenames = { "negative_z_e6",
                                           "negative_z_e7",
                                           "negative_z_e8",
                                           "negative_z_e9",
                                           "negative_z_e10",
                                           "negative_z_e11",
                                           "negative_z_e12" };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-1.5, -1.5, -1.5};
        PointType upper_bound = {1.5, 1.5, 1.5};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {0.0, 0.0, -perturbations[i]};
        RunCube(delta, lower_bound, upper_bound, perturbation, filenames[i]);
    }
}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra
