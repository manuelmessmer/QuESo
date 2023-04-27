// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// STL includes
#include <omp.h>
#include <string>
//// Project includes
#include "containers/triangle_mesh.hpp"
#include "io/io_utilities.h"
#include "embedding/brep_operator.h"
#include "utilities/mesh_utilities.h"
#include "tests/cpp_tests/class_testers/flood_fill_tester.hpp"

namespace tibra {
namespace Testing{


BOOST_AUTO_TEST_SUITE( Thingi10KTestSuite )

BOOST_AUTO_TEST_CASE( ElementClassificationTester ) {

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
    for( int i = 0; i < static_cast<int>(filenames.size()); ++i ){
        std::string filename = filenames[i];
        std::cout << filename << std::endl;
        // Read triangle mesh
        TriangleMesh triangle_mesh{};
        IO::ReadMeshFromSTL(triangle_mesh, filename.c_str());

        // Get min/max of triangle mesh
        auto bounding_box_stl = MeshUtilities::BoundingBox(triangle_mesh);
        PointType lower_bound_stl = bounding_box_stl.first;
        PointType upper_bound_stl = bounding_box_stl.second;
        std::cout << count << std::endl;
        count++;
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
        FloodFillTester filler(&brep_operator, parameters);
        auto result = filler.ClassifyElementsForTest();

        auto p_states = std::move(result.first);
        auto p_groups = std::move(result.second);

        std::pair<Unique<FloodFill::StatusVectorType>, Unique<FloodFill::GroupVectorSetType>> results_single;


        omp_set_num_threads(1);
        results_single = filler.ClassifyElementsForTest();
        omp_set_num_threads(omp_get_num_procs());
        auto p_states_single = std::move(results_single.first);
        auto p_groups_single = std::move(results_single.second);

        BOOST_CHECK_EQUAL( p_groups->size(), p_groups_single->size() );

        std::sort(p_groups->begin(), p_groups->end(), [](auto &rLHs, auto &rRHS) -> bool
            { if( std::get<1>(rLHs) == std::get<1>(rRHS) ) { return std::get<2>(rLHs) < std::get<2>(rRHS); }
              else { return std::get<1>(rLHs) < std::get<1>(rRHS); } });

        std::sort(p_groups_single->begin(), p_groups_single->end(), [](auto &rLHs, auto &rRHS) -> bool
            { if( std::get<1>(rLHs) == std::get<1>(rRHS) ) { return std::get<2>(rLHs) < std::get<2>(rRHS); }
            else { return std::get<1>(rLHs) < std::get<1>(rRHS); } });

        for( IndexType i = 0; i < p_groups->size(); ++i){
            auto group_size = std::get<1>((*p_groups)[i]).size();
            auto group_size_single = std::get<1>((*p_groups_single)[i]).size();
            BOOST_CHECK_EQUAL(group_size, group_size_single);
            auto group_inside_count = std::get<2>((*p_groups)[i]);
            auto group_inside_count_single = std::get<2>((*p_groups_single)[i]);
            BOOST_CHECK_EQUAL(group_inside_count, group_inside_count_single);
        }

        IndexType count_2 = 0;
        Mapper mapper(parameters);

        #pragma omp parallel for
        for(IndexType i = 0; i < num_elements[0]; ++i ){
            for(IndexType j = 0; j < num_elements[1]; ++j ){
                for(IndexType k = 0; k < num_elements[2]; ++k ){
                    IndexType index = mapper.GetVectorIndexFromMatrixIndices(i, j, k);
                    auto box = mapper.GetBoundingBoxFromIndex(i, j, k);
                    Vector3d local_lower_bound = box.first;
                    Vector3d local_upper_bound = box.second;

                    auto local_lower_bound_param = mapper.PointFromGlobalToParam(local_lower_bound);
                    auto local_upper_bound_param = mapper.PointFromGlobalToParam(local_upper_bound);

                    Element element(1, local_lower_bound_param, local_upper_bound_param, parameters);

                    auto p_clipped_mesh = brep_operator.pClipTriangleMeshUnique(local_lower_bound, local_upper_bound );
                    test_area += MeshUtilities::Area(*p_clipped_mesh);

                    auto status = brep_operator.GetIntersectionState(element);
                    if( status != (*p_states)[index] ){
                        std::cout << "count: " << count_2 << std::endl;
                        std::cout << "status: " << status << ", " << (*p_states)[index] << std::endl;
                        auto p_cube = MeshUtilities::pGetCuboid(local_lower_bound, local_upper_bound);
                        //auto p_clip = brep_operator.pClipTriangleMesh(local_lower_bound, local_upper_bound);
                       // IO::WriteMeshToSTL(*p_clip, "p_clip.stl", true);
                        IO::WriteMeshToSTL(*p_cube, "p_cube.stl", true);
                        TIBRA_ERROR("Here: ") << "Yes!! " << filename << std::endl;
                    }
                }
            }
        }
    }

    std::cout << "Successfull tests: " << count << std::endl;
}


BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra