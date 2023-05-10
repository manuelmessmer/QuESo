// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "containers/triangle_mesh.hpp"
#include "io/io_utilities.h"
#include "embedding/brep_operator.h"
#include "quadrature/trimmed_element.h"
#include "tests/cpp_tests/class_testers/trimmed_element_tester.hpp"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( TrimmedDomainTestSuite )

BOOST_AUTO_TEST_CASE(TrimemdDomainElephantTest) {
    TIBRA_INFO << "Testing :: Test Trimmed Domain :: Elephant" << std::endl;

    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};
    Vector3d lower_bound = {-0.4, -0.6, -0.35 };
    Vector3d upper_bound = {0.4, 0.6, 0.35 };
    Parameters parameters( {Component("lower_bound", lower_bound),
                            Component("upper_bound", upper_bound),
                            Component("min_num_boundary_triangles", 200UL),
                            Component("number_of_elements", number_of_elements),
                            Component("min_element_volume_ratio", 0.0),
                            Component("polynomial_order", Vector3i(2,2,2) ) } );



    auto p_triangle_mesh = TriangleMesh::New();
    IO::ReadMeshFromSTL(*p_triangle_mesh, "tibra/tests/cpp_tests/data/elephant.stl");

    // Build brep_operator
    BRepOperator brep_operator(p_triangle_mesh,parameters);
    const auto& r_triangle_mesh = brep_operator.GetTriangleMesh();

    const double delta_x = 0.1;
    const double delta_y = 0.1;
    const double delta_z = 0.1;

    std::ifstream file("tibra/tests/cpp_tests/results/surface_integral_elephant.txt");
    std::string line{};

    const double volume_ref = MeshUtilities::Volume(r_triangle_mesh);
    double volume_test = 0.0;
    const double area_ref = MeshUtilities::Area(r_triangle_mesh);
    double area_test = 0.0;

    IndexType number_trimmed_elements = 0;
    for(double x = lower_bound[0]; x <= upper_bound[0]; x += delta_x){
        for(double y = lower_bound[1]; y <= upper_bound[1]; y += delta_y){
            for(double z = lower_bound[2]; z <= upper_bound[2]; z += delta_z){
                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto local_lower_bound_param = Mapping::PointFromGlobalToParam(local_lower_bound, lower_bound, upper_bound);
                auto local_upper_bound_param = Mapping::PointFromGlobalToParam(local_upper_bound, lower_bound, upper_bound);

                auto p_clipped_mesh = brep_operator.pClipTriangleMeshUnique(local_lower_bound, local_upper_bound);
                area_test += MeshUtilities::Area(*p_clipped_mesh);

                Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);
                const auto status = brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound );
                if( status == IntersectionStatus::Trimmed){
                    // Get trimmed domain
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(local_lower_bound, local_upper_bound);

                    // Get volume
                    const auto& r_mesh = p_trimmed_domain->GetTriangleMesh();
                    volume_test += MeshUtilities::Volume(r_mesh);

                    // Get boundary integration points
                    auto p_boundary_ips = p_trimmed_domain->pGetBoundaryIps();

                    VectorType constant_terms{};
                    QuadratureTrimmedElementTester::ComputeConstantTerms(constant_terms, p_boundary_ips, element, parameters);

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
                    double norm_ref = 0.0;
                    for( auto& value : constant_terms ){
                        getline(file, line);
                        double ref_value = std::stod(line);
                        error += std::abs(value-ref_value);
                        norm_ref += std::abs(ref_value);
                        norm += std::abs(value);
                    }

                    if( norm_ref/constant_terms.size() > 1e-12 ){
                        BOOST_CHECK_LT( error/norm_ref, 1e-6);
                    } else {
                        BOOST_CHECK_LT( norm/constant_terms.size(), 1e-12);
                    }
                    number_trimmed_elements++;
                }
                else if( status == IntersectionStatus::Inside ){
                    volume_test += (local_upper_bound[0]-local_lower_bound[0])*(local_upper_bound[1]-local_lower_bound[1])*(local_upper_bound[2]-local_lower_bound[2]);
                }
            }
        }
    }
    file.close();
    BOOST_CHECK_LT( std::abs(area_test-area_ref)/area_ref, 1e-12);
    BOOST_CHECK_LT( std::abs(volume_test-volume_ref)/volume_ref, 1e-12);
    BOOST_CHECK_EQUAL(number_trimmed_elements, 166);

}

BOOST_AUTO_TEST_CASE(TrimmedDomainBunnyTest) {
    TIBRA_INFO << "Testing :: Test Trimmed Domain :: Bunny" << std::endl;

    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};
    Vector3d lower_bound = {-24.0, -43.0, 5.0 };
    Vector3d upper_bound = {85, 46.0, 115 };
    Parameters parameters( {Component("lower_bound", lower_bound),
                            Component("upper_bound", upper_bound),
                            Component("number_of_elements", number_of_elements),
                            Component("min_num_boundary_triangles", 100UL),
                            Component("min_element_volume_ratio", 0.0),
                            Component("polynomial_order", Vector3i(2,2,2) ) } );
    auto p_triangle_mesh = TriangleMesh::New();
    IO::ReadMeshFromSTL(*p_triangle_mesh, "tibra/tests/cpp_tests/data/stanford_bunny.stl");

    // Build brep_operator
    BRepOperator brep_operator(p_triangle_mesh, parameters);
    const auto& r_triangle_mesh = brep_operator.GetTriangleMesh();
    const double delta_x = 15;
    const double delta_y = 15;
    const double delta_z = 15;

    std::ifstream file("tibra/tests/cpp_tests/results/surface_integral_bunny.txt");
    std::string line{};

    const double volume_ref = MeshUtilities::Volume(r_triangle_mesh);
    double volume_test = 0.0;
    const double area_ref = MeshUtilities::Area(r_triangle_mesh);
    double area_test = 0.0;

    IndexType number_trimmed_elements = 0;
    for(double x = lower_bound[0]; x <= upper_bound[0]; x += delta_x){
        for(double y = lower_bound[1]; y <= upper_bound[1]; y += delta_y){
            for(double z = lower_bound[2]; z <= upper_bound[2]; z += delta_z){
                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto local_lower_bound_param = Mapping::PointFromGlobalToParam(local_lower_bound, lower_bound, upper_bound);
                auto local_upper_bound_param = Mapping::PointFromGlobalToParam(local_upper_bound, lower_bound, upper_bound);

                auto p_clipped_mesh = brep_operator.pClipTriangleMeshUnique(local_lower_bound, local_upper_bound);
                area_test += MeshUtilities::Area(*p_clipped_mesh);

                Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);

                const auto status = brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound );
                if( status == IntersectionStatus::Trimmed){
                    // Get Trimmed domain
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(local_lower_bound, local_upper_bound);

                    // Get volume
                    const auto& r_mesh = p_trimmed_domain->GetTriangleMesh();
                    volume_test += MeshUtilities::Volume(r_mesh);

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
                    QuadratureTrimmedElementTester::ComputeConstantTerms(constant_terms, p_boundary_ips, element, parameters);

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
                } else if( status == IntersectionStatus::Inside ){
                    volume_test += (local_upper_bound[0]-local_lower_bound[0])*(local_upper_bound[1]-local_lower_bound[1])*(local_upper_bound[2]-local_lower_bound[2]);
                }
            }
        }
    }
    file.close();

    BOOST_CHECK_LT( std::abs(area_test-area_ref)/area_ref, 1e-12);
    BOOST_CHECK_LT( std::abs(volume_test-volume_ref)/volume_ref, 1e-12);
    BOOST_CHECK_EQUAL(number_trimmed_elements, 171);
}

BOOST_AUTO_TEST_CASE(TestTrimmedDomainCylinderTest) {
    TIBRA_INFO << "Testing :: Test Trimmed Domain :: Cylinder" << std::endl;

    typedef boost::numeric::ublas::vector<double> VectorType;

    Vector3i number_of_elements = {1, 1, 1};

    Vector3d lower_bound = {-1.5, -1.5, -1 };
    Vector3d upper_bound = {1.5, 1.5, 12 };
    Parameters parameters( {Component("lower_bound", lower_bound),
                            Component("upper_bound", upper_bound),
                            Component("number_of_elements", number_of_elements),
                            Component("min_num_boundary_triangles", 100UL),
                            Component("min_element_volume_ratio", 0.0),
                            Component("polynomial_order", Vector3i(2,2,2) ) } );

    auto p_triangle_mesh = TriangleMesh::New();
    IO::ReadMeshFromSTL(*p_triangle_mesh, "tibra/tests/cpp_tests/data/cylinder.stl");

    // Build brep_operator
    BRepOperator brep_operator(p_triangle_mesh, parameters);
    const auto& r_triangle_mesh = brep_operator.GetTriangleMesh();

    const double delta_x = 1;
    const double delta_y = 1;
    const double delta_z = 1;

    std::ifstream file("tibra/tests/cpp_tests/results/surface_integral_cylinder.txt");
    std::string line{};

    const double volume_ref = MeshUtilities::Volume(r_triangle_mesh);
    double volume_test = 0.0;
    const double area_ref = MeshUtilities::Area(r_triangle_mesh);
    double area_test = 0.0;

    IndexType number_trimmed_elements = 0;
    for(double x = lower_bound[0]; x <= upper_bound[0]; x += delta_x){
        for(double y = lower_bound[1]; y <= upper_bound[1]; y += delta_y){
            for(double z = lower_bound[2]; z <= upper_bound[2]; z += delta_z){
                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto local_lower_bound_param = Mapping::PointFromGlobalToParam(local_lower_bound, lower_bound, upper_bound);
                auto local_upper_bound_param = Mapping::PointFromGlobalToParam(local_upper_bound, lower_bound, upper_bound);

                auto p_clipped_mesh = brep_operator.pClipTriangleMeshUnique(local_lower_bound, local_upper_bound);
                area_test += MeshUtilities::Area(*p_clipped_mesh);

                Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);

                const auto status = brep_operator.GetIntersectionState(local_lower_bound, local_upper_bound);
                if( status == IntersectionStatus::Trimmed){
                    // Get Trimmed domain
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(local_lower_bound, local_upper_bound);

                    // Get volume
                    const auto& r_mesh = p_trimmed_domain->GetTriangleMesh();
                    volume_test += MeshUtilities::Volume(r_mesh);

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
                    QuadratureTrimmedElementTester::ComputeConstantTerms(constant_terms, p_boundary_ips, element, parameters);

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
                } else if( status == IntersectionStatus::Inside ){
                    volume_test += (local_upper_bound[0]-local_lower_bound[0])*(local_upper_bound[1]-local_lower_bound[1])*(local_upper_bound[2]-local_lower_bound[2]);
                }
            }
        }
    }
    file.close();

    BOOST_CHECK_LT( std::abs(area_test-area_ref)/area_ref, 1e-12);
    BOOST_CHECK_LT( std::abs(volume_test-volume_ref)/volume_ref, 1e-12);
    BOOST_CHECK_EQUAL( number_trimmed_elements, 80);
}


void RunCubeWithCavity(const PointType rDelta, const PointType rLowerBound, const PointType rUpperBound,
    const PointType Perturbation ){

    Vector3i number_of_elements = {1, 1, 1};
    Parameters parameters( {Component("lower_bound", rLowerBound),
                            Component("upper_bound", rUpperBound),
                            Component("number_of_elements", number_of_elements),
                            Component("min_num_boundary_triangles", 100UL),
                            Component("min_element_volume_ratio", 0.0),
                            Component("polynomial_order", Vector3i(2,2,2) ) } );

    auto p_triangle_mesh = TriangleMesh::New();
    IO::ReadMeshFromSTL(*p_triangle_mesh, "tibra/tests/cpp_tests/data/cube_with_cavity.stl");

    auto& vertices = p_triangle_mesh->GetVertices();
    for( auto& v : vertices ){
        v[0] += Perturbation[0];
        v[1] += Perturbation[1];
        v[2] += Perturbation[2];
    }

    // Build brep_operator
    BRepOperator brep_operator(p_triangle_mesh, parameters);
    const auto& r_triangle_mesh = brep_operator.GetTriangleMesh();

    const double delta_x = rDelta[0];
    const double delta_y = rDelta[1];
    const double delta_z = rDelta[2];

    const double volume_ref = MeshUtilities::Volume(r_triangle_mesh);
    const double area_ref = MeshUtilities::Area(r_triangle_mesh);
    double volume = 0.0;
    double area = 0.0;
    IndexType number_trimmed_elements = 0;
    for(double x = rLowerBound[0]; x <= rUpperBound[0]; x += delta_x){
        for(double y = rLowerBound[1]; y <= rUpperBound[1]; y += delta_y){
            for(double z = rLowerBound[2]; z <= rUpperBound[2]; z += delta_z){
                Vector3d local_lower_bound = {x, y, z};
                Vector3d local_upper_bound = {x+delta_x, y+delta_y, z+delta_z};

                auto local_lower_bound_param = Mapping::PointFromGlobalToParam(local_lower_bound, rLowerBound, rUpperBound);
                auto local_upper_bound_param = Mapping::PointFromGlobalToParam(local_upper_bound, rLowerBound, rUpperBound);
                Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);

                auto p_clipped_mesh = brep_operator.pClipTriangleMeshUnique(local_lower_bound, local_upper_bound);
                area += MeshUtilities::Area(*p_clipped_mesh);
                // Get Trimmed domain
                auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(local_lower_bound, local_upper_bound);
                if( p_trimmed_domain ){
                    auto& r_mesh = p_trimmed_domain->GetTriangleMesh();
                    volume += MeshUtilities::Volume(r_mesh);
                    number_trimmed_elements++;
                }

            }
        }
    }

    BOOST_CHECK_LT( std::abs(area-area_ref)/area_ref, 1e-10 );
    BOOST_CHECK_LT( std::abs(volume-volume_ref)/volume_ref, 1e-9 );
}

BOOST_AUTO_TEST_CASE(TestTrimemdDomainCube1Test) {
    TIBRA_INFO << "Testing :: Test Trimmed Domain :: Cube 1" << std::endl;

    std::vector<double> perturbations = { 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16 };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-3.0, -3.0, -3.0};
        PointType upper_bound = {3.0, 3.0, 3.0};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {perturbations[i], 0.0 , 0.0};
        RunCubeWithCavity(delta, lower_bound, upper_bound, perturbation);
    }
}

BOOST_AUTO_TEST_CASE(TestTrimemdDomainCube2Test) {
    TIBRA_INFO << "Testing :: Test Trimmed Domain :: Cube 2" << std::endl;

    std::vector<double> perturbations = { 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16 };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-3.0, -3.0, -3.0};
        PointType upper_bound = {3.0, 3.0, 3.0};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {0.0, perturbations[i], 0.0};
        RunCubeWithCavity(delta, lower_bound, upper_bound, perturbation);
    }
}

BOOST_AUTO_TEST_CASE(TestTrimemdDomainCube3Test) {
    TIBRA_INFO << "Testing :: Test Trimmed Domain :: Cube 3" << std::endl;

    std::vector<double> perturbations = { 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16 };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-3.0, -3.0, -3.0};
        PointType upper_bound = {3.0, 3.0, 3.0};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {0.0, 0.0, -perturbations[i]};
        RunCubeWithCavity(delta, lower_bound, upper_bound, perturbation);
    }
}

BOOST_AUTO_TEST_CASE(TestTrimemdDomainCube4Test) {
    TIBRA_INFO << "Testing :: Test Trimmed Domain :: Cube 4" << std::endl;

    std::vector<double> perturbations = { 1e-6}; //, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16 };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-3.0, -3.0, -3.0};
        PointType upper_bound = {3.0, 3.0, 3.0};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {-perturbations[i], 0.0 , 0.0};
        RunCubeWithCavity(delta, lower_bound, upper_bound, perturbation);
    }
}

BOOST_AUTO_TEST_CASE(TestTrimemdDomainCube5Test) {
    TIBRA_INFO << "Testing :: Test Trimmed Domain :: Cube 5" << std::endl;

    std::vector<double> perturbations = { 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16 };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-3.0, -3.0, -3.0};
        PointType upper_bound = {3.0, 3.0, 3.0};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {0.0, -perturbations[i], 0.0};
        RunCubeWithCavity(delta, lower_bound, upper_bound, perturbation);
    }
}

BOOST_AUTO_TEST_CASE(TestTrimemdDomainCube6Test) {
    TIBRA_INFO << "Testing :: Test Trimmed Domain :: Cube 6" << std::endl;

    std::vector<double> perturbations = { 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16 };

    for( IndexType i = 0; i < perturbations.size(); ++i ){
        PointType lower_bound = {-3.0, -3.0, -3.0};
        PointType upper_bound = {3.0, 3.0, 3.0};
        PointType delta = {1.5, 1.5, 1.5};
        PointType perturbation = {0.0, 0.0, -perturbations[i]};
        RunCubeWithCavity(delta, lower_bound, upper_bound, perturbation);
    }
}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra
