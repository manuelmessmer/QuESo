// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "containers/triangle_mesh.hpp"
#include "quadrature/trimmed_element.h"
#include "tests/cpp_tests/class_testers/trimmed_element_tester.hpp"
#include "io/io_utilities.h"
#include "embedding/brep_operator.h"
#include "TIBRA_main.h"
#include "utilities/mesh_utilities.h"

namespace tibra {
namespace Testing{


BOOST_AUTO_TEST_SUITE( Thingi10KTestSuite )

BOOST_AUTO_TEST_CASE( Thingi10KFileTester ) {

    int argc_ = boost::unit_test::framework::master_test_suite().argc;
    char** argv_ = boost::unit_test::framework::master_test_suite().argv;

    TIBRA_ERROR_IF("Thingi10KFileTester", argc_ != 2 ) << "Please provide Filename STL: -- Filename \n";

    const std::string filename = argv_[1];

    const IndexType n_max = 50;
    const IndexType n_min = 5;
    int count = 0;
    // Read triangle mesh
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, filename.c_str());

    const double volume = MeshUtilities::Volume(triangle_mesh);

    // if( std::abs(volume - volume_2) / std::abs(volume) > 1 ){
    //     TIBRA_ERROR("here") << filename << std::endl;
    // }

    if( volume > EPS0 && MeshUtilities::EstimateQuality(triangle_mesh) <  1e-10 ){
        count++;
        // Get min/max of triangle mesh
        PointType lower_bound_stl = {MAXD, MAXD, MAXD};
        PointType upper_bound_stl = {LOWESTD, LOWESTD, LOWESTD};
        for( int i = 0; i < triangle_mesh.NumOfTriangles(); ++i){
            const auto& p1 = triangle_mesh.P1(i);
            const auto& p2 = triangle_mesh.P2(i);
            const auto& p3 = triangle_mesh.P3(i);

            const PointType x_values = {p1[0], p2[0], p3[0]};
            const PointType y_values = {p1[1], p2[1], p3[1]};
            const PointType z_values = {p1[2], p2[2], p3[2]};

            auto x_min_max = std::minmax_element(x_values.begin(), x_values.end());
            auto y_min_max = std::minmax_element(y_values.begin(), y_values.end());
            auto z_min_max = std::minmax_element(z_values.begin(), z_values.end());

            lower_bound_stl[0] = std::min<double>(*x_min_max.first, lower_bound_stl[0]);
            upper_bound_stl[0] = std::max<double>(*x_min_max.second, upper_bound_stl[0]);

            lower_bound_stl[1] = std::min<double>(*y_min_max.first, lower_bound_stl[1]);
            upper_bound_stl[1] = std::max<double>(*y_min_max.second, upper_bound_stl[1]);

            lower_bound_stl[2] = std::min<double>(*z_min_max.first, lower_bound_stl[2]);
            upper_bound_stl[2] = std::max<double>(*z_min_max.second, upper_bound_stl[2]);
        }


        auto delta_stl = upper_bound_stl - lower_bound_stl;
        double h = 1.2*std::min( Math::Max(delta_stl)/n_max, Math::Min(delta_stl)/n_min );

        PointType lower_bound = lower_bound_stl - delta_stl*0.1;

        Vector3i num_elements{};
        num_elements[0] = std::ceil( 1.2* delta_stl[0] / h );
        num_elements[1] = std::ceil( 1.2* delta_stl[1] / h );
        num_elements[2] = std::ceil( 1.2* delta_stl[2] / h );

        PointType upper_bound{};
        upper_bound[0] = lower_bound[0] + num_elements[0]*h;
        upper_bound[1] = lower_bound[1] + num_elements[1]*h;
        upper_bound[2] = lower_bound[2] + num_elements[2]*h;

        Parameters parameters( {
                    Component("input_filename", filename),
                    Component("lower_bound", lower_bound),
                    Component("echo_level", 3UL),
                    Component("upper_bound", upper_bound),
                    Component("min_num_boundary_triangles", 10UL),
                    Component("number_of_elements", num_elements),
                    Component("min_element_volume_ratio", 0.0),
                    Component("integration_method", IntegrationMethod::Gauss),
                    Component("moment_fitting_residual", 1e-6),
                    Component("polynomial_order", Vector3i(2,2,2) ) } );

        auto delta_new = upper_bound - lower_bound;

        double test_volume = 0.0;
        double test_volume_ips = 0.0;
        double test_area = 0.0;
        double vol_inner = 0.0;
        std::vector<double> volumes;
        BRepOperator brep_operator(triangle_mesh, parameters);
        IndexType id_el = 0;

        IndexType number_trimmed_el = 0;
        TriangleMesh clipped_meshes;
        clipped_meshes.Reserve(5);
        for(double x = lower_bound[0]; x <= upper_bound[0]-EPS0; x += h){
            for(double y = lower_bound[1]; y <= upper_bound[1]-EPS0; y += h){
                for(double z = lower_bound[2]; z <= upper_bound[2]-EPS0; z += h){
                    Vector3d local_lower_bound = {x, y, z};
                    Vector3d local_upper_bound = {x+h, y+h, z+h};

                    auto local_lower_bound_param = Mapping::GlobalToParam(local_lower_bound, lower_bound, upper_bound);
                    auto local_upper_bound_param = Mapping::GlobalToParam(local_upper_bound, lower_bound, upper_bound);

                    Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);

                    auto p_clipped_mesh = brep_operator.pClipTriangleMeshUnique(local_lower_bound, local_upper_bound );
                    test_area += MeshUtilities::Area(*p_clipped_mesh);

                    auto status = brep_operator.GetIntersectionState(element);
                    if( status == IntersectionStatus::Trimmed){
                        number_trimmed_el++;
                        // Get trimmed domain
                        auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(local_lower_bound, local_upper_bound);
                        if( p_trimmed_domain ){
                            // Add trimmed domain to element.
                            auto mesh = p_trimmed_domain->GetTriangleMesh();
                            test_volume += MeshUtilities::Volume(mesh);
                            // element.SetIsTrimmed(true);
                            // element.pSetTrimmedDomain(p_trimmed_domain);
                            // const auto residual = QuadratureTrimmedElementTester::AssembleIPs(element, parameters);
                            // auto& r_points = element.GetIntegrationPoints();
                            // for( auto point : r_points ){
                            //     test_volume_ips += point.GetWeight()*delta_new[0]*delta_new[1]*delta_new[2];
                            // }
                        }
                    } else if( status == IntersectionStatus::Inside ){
                        test_volume += (local_upper_bound[0]-local_lower_bound[0]) * (local_upper_bound[1]-
                        local_lower_bound[1]) * (local_upper_bound[2] - local_lower_bound[2]);
                        // test_volume_ips += (local_upper_bound[0]-local_lower_bound[0]) * (local_upper_bound[1]-
                        // local_lower_bound[1]) * (local_upper_bound[2] - local_lower_bound[2]);

                    }
                }
            }
        }

        const double area = MeshUtilities::Area(triangle_mesh);
        const double error_area = std::abs(area - test_area)/area;
        const double error_volume = std::abs( test_volume - volume)/ volume;

    }

}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra