// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef PARAMETERS_INCLUDE_H
#define PARAMETERS_INCLUDE_H

//// STL includes
#include <variant>
#include <iostream>
#include <string>
#include <exception>
#include <type_traits>
#include <vector>
////Project includes
#include "containers/point_types.h"
//#include "quadrature/integration_points_1d/integration_points_factory_1d.h"


namespace tibra {

///@}
///@name  Tibra Globals
///@{

enum IntegrationMethod {Gauss, ReducedGauss1, ReducedGauss2, ReducedExact, ReducedOrder1, ReducedOrder2};
typedef enum IntegrationMethod IntegrationMethodType;


class Component {
public:
    //typedef IntegrationPointFactory1D::IntegrationMethod IntegrationMethod;
    typedef std::variant<PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType> ComponentType;

    // Constructor
    Component(std::string Name, ComponentType Component) : mName(Name), mComponents(Component)
    {}

    const ComponentType& Get() const{
        return mComponents;
    }

    const std::string& Name() const {
        return mName;
    }

private:
    std::string mName{};
    ComponentType mComponents{};
};


class Parameters {
public:
    //typedef IntegrationPointFactory1D::IntegrationMethod IntegrationMethod;
    struct TypeVisit {
        const std::type_info* operator()(const PointType& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const Vector3i& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const IndexType& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const double& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const std::string& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const bool& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const IntegrationMethodType& rValue){return &typeid(rValue); };
    };

    struct PrintVisit {
        PrintVisit(std::ostream& rOStream) : mOStream(rOStream){}
        void operator()(const PointType& rValue){mOStream << rValue; };
        void operator()(const Vector3i& rValue){mOStream << rValue;};
        void operator()(const IndexType& rValue){ mOStream << rValue;};
        void operator()(const double& rValue){mOStream << rValue;};
        void operator()(const std::string& rValue){mOStream << rValue;};
        void operator()(const bool& rValue){mOStream << rValue; };
        void operator()(const IntegrationMethodType& rValue){ mOStream << rValue; };

    private:
        std::ostream& mOStream;
    };

    Parameters() {
        AddDefaults();
        CheckTypes();
    }

    Parameters(std::vector<Component> Component) : mComponents(Component)
    {
        AddDefaults();
        CheckTypes();
    }

    const PointType& LowerBound() const {
        return Get<PointType>("lower_bound");
    }

    const PointType& UpperBound() const {
        return Get<PointType>("lower_bound");
    }

    const PointType& Order() const {
        return Get<PointType>("polynomial_order");
    }

    const IntegrationMethodType& IntegrationMethod() const {
        return Get<IntegrationMethodType>("integration_method");
    }

    const Vector3i& NumberOfElements() const {
        return Get<Vector3i>("number_of_knot_spans");;
    }

    double InitialTriangleEdgeLength() const {
        /// deprecated
        return Get<double>("initial_triangle_edge_length");
    }

    IndexType MinimumNumberOfTriangles() const {
        /// deprecated
        return Get<IndexType>("min_num_boundary_triangles");
    }

    double MomentFittingResidual() const {
        return Get<double>("moment_fitting_residual");
    }

    IndexType EchoLevel() const {
        return Get<IndexType>("echo_level");
    }

    bool UseCustomizedTrimmedPositions() const{
        /// Only used for testing.
        return Get<bool>("use_customized_trimmed_points");
    }

    IndexType GetPointDistributionFactor() const {
        return Get<IndexType>("init_point_distribution_factor");
    }

    bool GGQRuleIsUsed() const {
        return IntegrationMethod() >= 3;
    }

    template<typename type>
    void Set(const std::string& rName, const type& rValue){
        auto p_value = pFind<type>(rName);
        if( p_value ){
            *p_value = rValue;
        }

        mComponents.push_back(Component(rName, rValue));
        CheckTypes();
    }

    template<typename type>
    const type& Get( std::string Name ) const {
        const auto p_value = pFind<type>(Name);
        if( p_value ){
            return *p_value;
        }

        std::string error_message = "Parameter :: Get :: Component: '" + Name + "' not found.\n";
        throw std::runtime_error(error_message);
    }

    void PrintInfo(std::ostream& rOStream) const {
        rOStream << "Parameters: \n";
        for( auto& value : mComponents ){
            rOStream << value.Name() << ": ";
            std::visit(PrintVisit(rOStream), value.Get());
            rOStream << "\n";
        }
    }

private:


    void AddDefaults(){
        for( auto& value : mDefaults){
            AddValueNoThrow<PointType>(value);
            AddValueNoThrow<Vector3i>(value);
            AddValueNoThrow<bool>(value);
            AddValueNoThrow<double>(value);
            AddValueNoThrow<IndexType>(value);
            AddValueNoThrow<std::string>(value);
            AddValueNoThrow<IntegrationMethodType>(value);
        }
    }

    template<typename type>
    void AddValueNoThrow(const Component& rValues){
        auto p_type_id = std::visit(TypeVisit{}, rValues.Get());
        if( *p_type_id == typeid(type) ){
            auto p_value = pFind<type>(rValues.Name());
            if( !p_value ){
                mComponents.push_back( rValues );
            }
        }
    }

    void CheckTypes() const {
        for( auto& r_components : mComponents ){
            const auto& current_name = r_components.Name();
            const auto& current_value = r_components.Get();
            IndexType count = 0;
            for( auto& str_type_pair : mTypes){
                const auto& ref_name = str_type_pair.first;
                const auto& p_ref_type_info = str_type_pair.second;
                if( ref_name ==  current_name ){
                    auto p_curren_type_info = std::visit(TypeVisit{}, r_components.Get());
                    if( *p_ref_type_info ==  *p_curren_type_info){
                        break;
                    } else {
                        std::string error_message = "Parameters :: CheckTypes :: Name: '" + current_name +
                            "' is not provided with correct Type.\n";
                        throw std::runtime_error(error_message);
                    }
                }
                count++;
            }
            if( count >= mTypes.size() ){
                std::string error_message = "Parameters :: CheckTypes :: Name: '" + current_name +
                    "' is not a valid Parameter.\n";
                throw std::runtime_error(error_message);
            }
        }
    }

    template<class type>
    const type* pFind( const std::string& rName ) const {
        for( auto& r_component : mComponents){
            const auto p_value = std::get_if<type>(&r_component.Get());
            if( p_value ){
                if( r_component.Name() == rName ){
                    return p_value;
                }
            }
        }
        return nullptr;
    }

    template<class type>
    type* pFind( const std::string& rName ) {
        for( auto& r_component : mComponents){
            auto p_value = const_cast<type*>(std::get_if<type>(&r_component.Get()));
            if( p_value ){
                if( r_component.Name() == rName ){
                    return p_value;
                }
            }
        }
        return nullptr;
    }

    std::vector<Component> mComponents{};
    std::vector<Component> mDefaults = { Component("echo_level", 0UL),
                                         Component("embedding_flag", true),
                                         Component("initial_triangle_edge_length", 1.0),
                                         Component("min_num_boundary_triangles", 1000UL),
                                         Component("moment_fitting_residual", 1.0e-10),
                                         Component("init_point_distribution_factor", 2UL),
                                         Component("polynomial_order", Vector3i(2UL, 2Ul, 2UL) ),
                                         Component("integration_method", IntegrationMethod::Gauss),
                                         Component("use_customized_trimmed_points", false) };

    std::vector<std::pair<std::string, const std::type_info*>> mTypes = {
       std::make_pair<std::string, const std::type_info*>("input_filename", &typeid(std::string) ),
       std::make_pair<std::string, const std::type_info*>("postprocess_filename", &typeid(std::string) ),
       std::make_pair<std::string, const std::type_info*>("echo_level", &typeid(IndexType) ),
       std::make_pair<std::string, const std::type_info*>("embedding_flag", &typeid(bool) ),
       std::make_pair<std::string, const std::type_info*>("lower_bound", &typeid(PointType) ),
       std::make_pair<std::string, const std::type_info*>("upper_bound", &typeid(PointType) ),
       std::make_pair<std::string, const std::type_info*>("polynomial_order", &typeid(Vector3i) ),
       std::make_pair<std::string, const std::type_info*>("number_of_knot_spans", &typeid(Vector3i) ),
       std::make_pair<std::string, const std::type_info*>("initial_triangle_edge_length", &typeid(double) ),
       std::make_pair<std::string, const std::type_info*>("min_num_boundary_triangles", &typeid(IndexType) ),
       std::make_pair<std::string, const std::type_info*>("moment_fitting_residual", &typeid(double) ),
       std::make_pair<std::string, const std::type_info*>("init_point_distribution_factor", &typeid(IndexType) ),
       std::make_pair<std::string, const std::type_info*>("integration_method", &typeid(IntegrationMethodType) ),
       std::make_pair<std::string, const std::type_info*>("use_customized_trimmed_points", &typeid(bool) ) };

}; // End class Parameters

std::ostream& operator<< (std::ostream& rOStream, const Parameters& rThis);

} // End namespace tibra

#endif // PARAMETERS_INCLUDE_H
