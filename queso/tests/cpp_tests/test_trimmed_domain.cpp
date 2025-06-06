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
#include "queso/includes/dictionary_factory.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/io/io_utilities.h"
#include "queso/embedding/brep_operator.h"
#include "queso/quadrature/trimmed_element.hpp"
#include "queso/tests/cpp_tests/class_testers/trimmed_element_tester.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( TrimmedDomainTestSuite )

void CheckTriangleOrientation(const TriangleMeshInterface& rTriangleMesh){
    for( IndexType i = 0; i < rTriangleMesh.NumOfTriangles(); ++i){
        if( rTriangleMesh.Area(i) > EPS3 ){
            const auto& p0 = rTriangleMesh.P1(i);
            const auto& p1 = rTriangleMesh.P2(i);
            const auto& p2 = rTriangleMesh.P3(i);
            const auto A = Math::Subtract(p1, p0);
            const auto B = Math::Subtract(p2, p1);

            PointType normal_cross = Math::Cross(A, B);
            Math::DivideSelf(normal_cross, Math::Norm(normal_cross) );

            PointType normal_stored = rTriangleMesh.Normal(i);

            QuESo_CHECK_POINT_NEAR(normal_cross, normal_stored, EPS2);
        }
    }
}

void RunTest(const std::string& rFilename, const Dictionary<queso::key::MainValuesTypeTag>& rSettings,
            const std::string& rResultsFilename, IndexType NumTrimmedElements){

    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, rFilename.c_str());

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    //std::ofstream file_out{};
    std::ifstream file(rResultsFilename);
    std::string line{};

    //file_out.open("test.txt");
    const double volume_ref = MeshUtilities::Volume(triangle_mesh);
    double volume_test = 0.0;
    const double area_ref = MeshUtilities::Area(triangle_mesh);
    double area_test = 0.0;

    const double min_vol_ratio = rSettings[MainSettings::trimmed_quadrature_rule_settings].
        GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio);
    const IndexType min_num_triangles = rSettings[MainSettings::trimmed_quadrature_rule_settings].
        GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles);
    const Vector3i order =  rSettings[MainSettings::background_grid_settings].
        GetValue<Vector3i>(BackgroundGridSettings::polynomial_order);

    GridIndexer grid_indexer(rSettings);
    IndexType number_trimmed_elements = 0;
    for( IndexType i = 0; i < grid_indexer.NumberOfElements(); ++i){
        const BoundingBoxType bounding_box = grid_indexer.GetBoundingBoxXYZFromIndex(i);
        const Vector3d lower_bound_xyz = bounding_box.first;
        const Vector3d upper_bound_xyz = bounding_box.second;
        const Vector3d delta_xyz = Math::Subtract( upper_bound_xyz, lower_bound_xyz );
        const BoundingBoxType bounding_box_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(i);
        const Vector3d lower_bound_uvw = bounding_box_uvw.first;
        const Vector3d upper_bound_uvw = bounding_box_uvw.second;

        // Construct element
        ElementType element(1, MakeBox(lower_bound_xyz, upper_bound_xyz),
                               MakeBox(lower_bound_uvw, upper_bound_uvw));

        auto p_clipped_mesh = brep_operator.pClipTriangleMeshUnique(lower_bound_xyz, upper_bound_xyz);
        area_test += MeshUtilities::Area(*p_clipped_mesh);


        const auto status = brep_operator.GetIntersectionState(lower_bound_xyz, upper_bound_xyz );
        if( status == IntersectionState::trimmed){
            // Get trimmed domain
            auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound_xyz, upper_bound_xyz, min_vol_ratio, min_num_triangles);

            // Get triangle mesh
            const auto& r_mesh = p_trimmed_domain->GetTriangleMesh();

            // Check orientation
            CheckTriangleOrientation(r_mesh);

            // Get volume
            volume_test += MeshUtilities::Volume(r_mesh);

            // Get boundary integration points
            auto p_boundary_ips = p_trimmed_domain->pGetBoundaryIps<BoundaryIntegrationPointType>();

            std::vector<double> constant_terms{};
            QuadratureTrimmedElementTester<ElementType>::ComputeConstantTerms(constant_terms, p_boundary_ips, element, order);

            // Read and ignore header
            getline(file, line);
            //file_out << "E: " << number_trimmed_elements << std::endl;
            double surface_area = 0.0;
            for( auto& ip : *p_boundary_ips){
                surface_area += ip.Weight();
            }
            //file_out << std::setprecision(16) << surface_area << std::endl;
            // Read ref surface area
            getline(file, line);
            double ref_surface_area = std::stod(line);
            QuESo_CHECK_RELATIVE_NEAR(surface_area, ref_surface_area, 1e-10);

            double error = 0.0;
            double norm = 0.0;
            double norm_ref = 0.0;
            for( auto& value : constant_terms ){
                //file_out << std::setprecision(16) << value << std::endl;
                getline(file, line);
                double ref_value = std::stod(line);
                error += std::abs(value-ref_value);
                norm_ref += std::abs(ref_value);
                norm += std::abs(value);
            }

            if( norm_ref/constant_terms.size() > 1e-12 ){
                QuESo_CHECK_LT( error/norm_ref, 1e-6);
            } else {
                QuESo_CHECK_LT( norm/constant_terms.size(), 1e-12);
            }
            number_trimmed_elements++;
        }
        else if( status == IntersectionState::inside ){
            volume_test += delta_xyz[0]*delta_xyz[1]*delta_xyz[2];
        }
    }
    //file_out.close();
    file.close();
    QuESo_CHECK_RELATIVE_NEAR(area_test, area_ref, 1e-12);
    QuESo_CHECK_RELATIVE_NEAR(volume_test, volume_ref, 1e-12);
    QuESo_CHECK_EQUAL(number_trimmed_elements, NumTrimmedElements);
}

BOOST_AUTO_TEST_CASE(TrimemdDomainElephantTest) {
    QuESo_INFO << "Testing :: Test Trimmed Domain :: Elephant" << std::endl;

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-0.4, -0.6, -0.35});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{0.5, 0.7, 0.45});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{0.0, 0.0, 0.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 1.0, 1.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, Vector3i{9, 13, 8});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::polynomial_order, Vector3i{2,2,2});

    r_settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, 200u);
    r_settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio, 0.0);

    RunTest("queso/tests/cpp_tests/data/elephant.stl", r_settings,
            "queso/tests/cpp_tests/results/surface_integral_elephant.txt", 166);
}

BOOST_AUTO_TEST_CASE(TrimmedDomainBunnyTest) {
    QuESo_INFO << "Testing :: Test Trimmed Domain :: Bunny" << std::endl;

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-24.0, -43.0, 5.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{85, 46.0, 115});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{-24.0, -43.0, 5.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{85, 46.0, 115});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, Vector3i{8, 6, 8});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::polynomial_order, Vector3i{2,2,2});

    r_settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, 100u);
    r_settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio, 0.0);

    RunTest("queso/tests/cpp_tests/data/stanford_bunny.stl", r_settings,
            "queso/tests/cpp_tests/results/surface_integral_bunny.txt", 186);
}

BOOST_AUTO_TEST_CASE(TestTrimmedDomainCylinderTest) {
    QuESo_INFO << "Testing :: Test Trimmed Domain :: Cylinder" << std::endl;

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-1.5, -1.5, -1});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{1.5, 1.5, 12});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{-1.0, -1.0, -1.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 1.0, 1.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, Vector3i{3, 3, 13});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::polynomial_order, Vector3i{2,2,2});

    r_settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, 100u);
    r_settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio, 0.0);

    RunTest("queso/tests/cpp_tests/data/cylinder.stl", r_settings,
            "queso/tests/cpp_tests/results/surface_integral_cylinder.txt", 80);
}


void RunCubeWithCavity(const PointType rDelta, const PointType rLowerBound, const PointType rUpperBound,
    const PointType Perturbation ){

    Vector3i number_of_elements = {1, 1, 1};

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, rLowerBound);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, rUpperBound);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, rLowerBound);
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, rUpperBound);
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, number_of_elements);
    r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{2,2,2});

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cube_with_cavity.stl");

    auto& vertices = triangle_mesh.GetVertices();
    for( auto& v : vertices ){
        v[0] += Perturbation[0];
        v[1] += Perturbation[1];
        v[2] += Perturbation[2];
    }

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 0.0;
    const IndexType min_num_triangles = 100;

    const double volume_ref = MeshUtilities::Volume(triangle_mesh);
    const double area_ref = MeshUtilities::Area(triangle_mesh);
    double volume = 0.0;
    double area = 0.0;

    GridIndexer grid_indexer(r_settings);
    IndexType number_trimmed_elements = 0;
    for( IndexType i = 0; i < grid_indexer.NumberOfElements(); ++i){
        const BoundingBoxType bounding_box = grid_indexer.GetBoundingBoxXYZFromIndex(i);
        const Vector3d lower_bound_xyz = bounding_box.first;
        const Vector3d upper_bound_xyz = bounding_box.second;

        auto p_clipped_mesh = brep_operator.pClipTriangleMeshUnique(lower_bound_xyz, upper_bound_xyz);
        area += MeshUtilities::Area(*p_clipped_mesh);
        // Get Trimmed domain
        auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound_xyz, upper_bound_xyz, min_vol_ratio, min_num_triangles);
        if( p_trimmed_domain ){
            auto& r_mesh = p_trimmed_domain->GetTriangleMesh();
            // Check triangle orientations.
            CheckTriangleOrientation(r_mesh);

            volume += MeshUtilities::Volume(r_mesh);
            number_trimmed_elements++;
        }
    }

    QuESo_CHECK_RELATIVE_NEAR(area, area_ref, 1e-10 );
    QuESo_CHECK_RELATIVE_NEAR(volume, volume_ref, 1e-9 );
}

BOOST_AUTO_TEST_CASE(TestTrimemdDomainCube1Test) {
    QuESo_INFO << "Testing :: Test Trimmed Domain :: Cube 1" << std::endl;

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
    QuESo_INFO << "Testing :: Test Trimmed Domain :: Cube 2" << std::endl;

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
    QuESo_INFO << "Testing :: Test Trimmed Domain :: Cube 3" << std::endl;

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
    QuESo_INFO << "Testing :: Test Trimmed Domain :: Cube 4" << std::endl;

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
    QuESo_INFO << "Testing :: Test Trimmed Domain :: Cube 5" << std::endl;

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
    QuESo_INFO << "Testing :: Test Trimmed Domain :: Cube 6" << std::endl;

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
} // End namespace queso
