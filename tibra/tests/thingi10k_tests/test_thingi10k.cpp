// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// STL includes
#include <string>
//// Project includes
#include "containers/triangle_mesh.hpp"
#include "io/io_utilities.h"
#include "embedding/brep_operator.h"
#include "utilities/mesh_utilities.h"

namespace tibra {
namespace Testing{


BOOST_AUTO_TEST_SUITE( Thingi10KTestSuite )

BOOST_AUTO_TEST_CASE( Thingi10KTester ) {

    int argc_ = boost::unit_test::framework::master_test_suite().argc;
    char** argv_ = boost::unit_test::framework::master_test_suite().argv;

    TIBRA_ERROR_IF("Thingi10KFileTester", argc_ != 5 ) << "Please provide following arguments: -- single/all Filename/Directory n_min n_max\n";

    std::vector<std::string> filenames;
    filenames.reserve(5000);
    std::string option = argv_[1];
    if( option == "single"){
        filenames.push_back(argv_[2]);
    } else if( option == "all" ) {
        std::ifstream file("tibra/tests/thingi10k_tests/model_ids.txt");
        std::string buffer;
        while( std::getline(file, buffer) ){
            std::string name = argv_[2] +  buffer + ".stl";
            filenames.push_back( name );
        }
    } else {
        TIBRA_ERROR("Thingi10KTester") << "First argument must be single/all. Provided: " << argv_[1];
    }

    const IndexType n_min = static_cast<IndexType>(std::stoi(argv_[3]));
    const IndexType n_max = static_cast<IndexType>(std::stoi(argv_[4]));

    std::cout << "Testing " << filenames.size() << " STLs\n";
    IndexType count = 0;
    #pragma omp parallel for reduction(+ : count) schedule(dynamic)
    for( int i = 0; i < static_cast<int>(filenames.size()); ++i ){
        std::string filename = filenames[i];
        // Read triangle mesh
        TriangleMesh triangle_mesh{};
        IO::ReadMeshFromSTL(triangle_mesh, filename.c_str());

        // Get min/max of triangle mesh
        auto bounding_box_stl = MeshUtilities::BoundingBox(triangle_mesh);
        PointType lower_bound_stl = bounding_box_stl.first;
        PointType upper_bound_stl = bounding_box_stl.second;

        auto delta_stl = upper_bound_stl - lower_bound_stl;
        double h = 1.2*std::min( Math::Max(delta_stl)/n_max, Math::Min(delta_stl)/n_min );

        PointType lower_bound = lower_bound_stl - delta_stl*0.1;

        Vector3i num_elements{};
        num_elements[0] = static_cast<IndexType>(std::ceil( 1.2* delta_stl[0] / h ));
        num_elements[1] = static_cast<IndexType>(std::ceil( 1.2* delta_stl[1] / h ));
        num_elements[2] = static_cast<IndexType>(std::ceil( 1.2* delta_stl[2] / h ));

        PointType upper_bound{};
        upper_bound[0] = lower_bound[0] + num_elements[0]*h;
        upper_bound[1] = lower_bound[1] + num_elements[1]*h;
        upper_bound[2] = lower_bound[2] + num_elements[2]*h;

        Parameters parameters( { Component("input_filename", filename),
                                Component("lower_bound", lower_bound),
                                Component("upper_bound", upper_bound),
                                Component("min_num_boundary_triangles", 10UL),
                                Component("number_of_elements", num_elements),
                                Component("min_element_volume_ratio", 0.0) } );

        auto delta_new = upper_bound - lower_bound;

        double test_volume = 0.0;
        double test_area = 0.0;

        BRepOperator brep_operator(triangle_mesh, parameters);
        for(double x = lower_bound[0]; x < upper_bound[0]-EPS0; x += h){
            for(double y = lower_bound[1]; y < upper_bound[1]-EPS0; y += h){
                for(double z = lower_bound[2]; z < upper_bound[2]-EPS0; z += h){
                    Vector3d local_lower_bound = {x, y, z};
                    Vector3d local_upper_bound = {x+h, y+h, z+h};

                    auto local_lower_bound_param = Mapping::GlobalToParam(local_lower_bound, lower_bound, upper_bound);
                    auto local_upper_bound_param = Mapping::GlobalToParam(local_upper_bound, lower_bound, upper_bound);
                    Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);

                    auto p_clipped_mesh = brep_operator.pClipTriangleMeshUnique(local_lower_bound, local_upper_bound );
                    test_area += MeshUtilities::Area(*p_clipped_mesh);

                    auto status = brep_operator.GetIntersectionState(element);
                    if( status == IntersectionStatus::Trimmed){
                        // Get trimmed domain
                        auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(local_lower_bound, local_upper_bound);
                        if( p_trimmed_domain ){
                            auto mesh = p_trimmed_domain->GetTriangleMesh();
                            test_volume += MeshUtilities::Volume(mesh);
                        }
                    } else if( status == IntersectionStatus::Inside ){
                        test_volume += (local_upper_bound[0] - local_lower_bound[0])
                                    * (local_upper_bound[1] - local_lower_bound[1])
                                    * (local_upper_bound[2] - local_lower_bound[2]);
                    }
                }
            }
        }

        const double area_ref = MeshUtilities::Area(triangle_mesh);
        const double error_area = std::abs(test_area - area_ref)/area_ref;

        const double volume_ref = MeshUtilities::Volume(triangle_mesh);
        const double error_volume = std::abs( test_volume - volume_ref)/ volume_ref;

        BOOST_CHECK_LT(error_volume, 3e-10);
        BOOST_CHECK_LT(error_area, 1e-10);

        if( error_volume < 3e-10 && error_area < 1e-10 ){
            count++;
        }
    }

    std::cout << "Successfull tests: " << count << std::endl;
}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra