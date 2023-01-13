#include "utilities/parameters.h"

namespace tibra {

const std::vector<Component> Parameters::mDefaultComponents = {
    Component("echo_level", 0UL),
    Component("embedding_flag", true),
    Component("initial_triangle_edge_length", 1.0),
    Component("min_num_boundary_triangles", 500UL),
    Component("moment_fitting_residual", 1.0e-10),
    Component("min_element_volume_ratio", 1.0e-3),
    Component("init_point_distribution_factor", 2UL),
    Component("polynomial_order", Vector3i(2UL, 2Ul, 2UL) ),
    Component("integration_method", IntegrationMethod::GGQ_Optimal),
    Component("use_customized_trimmed_points", false) };

const std::vector<std::pair<std::string, const std::type_info*>> Parameters::mAllAvailableComponents = {
    std::make_pair<std::string, const std::type_info*>("input_filename", &typeid(std::string) ),
    std::make_pair<std::string, const std::type_info*>("postprocess_filename", &typeid(std::string) ),
    std::make_pair<std::string, const std::type_info*>("echo_level", &typeid(IndexType) ),
    std::make_pair<std::string, const std::type_info*>("embedding_flag", &typeid(bool) ),
    std::make_pair<std::string, const std::type_info*>("lower_bound", &typeid(PointType) ),
    std::make_pair<std::string, const std::type_info*>("upper_bound", &typeid(PointType) ),
    std::make_pair<std::string, const std::type_info*>("polynomial_order", &typeid(Vector3i) ),
    std::make_pair<std::string, const std::type_info*>("number_of_elements", &typeid(Vector3i) ),
    std::make_pair<std::string, const std::type_info*>("initial_triangle_edge_length", &typeid(double) ),
    std::make_pair<std::string, const std::type_info*>("min_num_boundary_triangles", &typeid(IndexType) ),
    std::make_pair<std::string, const std::type_info*>("moment_fitting_residual", &typeid(double) ),
    std::make_pair<std::string, const std::type_info*>("min_element_volume_ratio", &typeid(double) ),
    std::make_pair<std::string, const std::type_info*>("init_point_distribution_factor", &typeid(IndexType) ),
    std::make_pair<std::string, const std::type_info*>("integration_method", &typeid(IntegrationMethodType) ),
    std::make_pair<std::string, const std::type_info*>("use_customized_trimmed_points", &typeid(bool) ) };

std::ostream& operator<< (std::ostream& rOStream, const Parameters& rThis){
    rThis.PrintInfo(rOStream);
    return rOStream;
}

} // End tibra namespace