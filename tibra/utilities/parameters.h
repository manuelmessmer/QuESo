// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef PARAMETERS_INCLUDE_H
#define PARAMETERS_INCLUDE_H

//// STL includes
#include <variant>
#include <iostream>
#include <string>
#include <algorithm>
#include <exception>
#include <type_traits>
#include <vector>
////Project includes
#include "define.hpp"

namespace tibra {

///@name  TIBRA Globals
///@{

enum IntegrationMethod {Gauss, Gauss_Reduced1, Gauss_Reduced2, GGQ_Optimal, GGQ_Reduced1, GGQ_Reduced2};
typedef enum IntegrationMethod IntegrationMethodType;

///@}
///@name  TIBRA Classes
///@{

/**
 * @class  Component
 * @author Manuel Messmer
 * @brief  Container for all available Types of parameters. Available Types are:
 *         PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType
 * @details Stores Name (Key) and Value of Parameter.
**/
class Component {
public:
    ///@name Type Definitions
    ///@{
    typedef std::variant<PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType> ComponentType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    Component(std::string Name, ComponentType NewComponent) : mName(Name), mComponents(NewComponent)
    {}

    ///@}
    ///@name Operations
    ///@{

    /// Returns Value
    const ComponentType& Get() const{
        return mComponents;
    }

    /// Returns Name (Key)
    const std::string& Name() const {
        return mName;
    }


private:
    ///@}
    ///@name Private Member variables
    ///@{

    std::string mName{};
    ComponentType mComponents{};

    ///@}

}; // End class Component


/**
 * @class  Parameters
 * @author Manuel Messmer
 * @brief  Dynamic container for all available parameters.
**/
class Parameters {
public:

    ///@name std::visit Structs
    ///@{
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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    Parameters() {
        AddDefaults();
        CheckComponents();
    }

    /// Constructor
    Parameters(std::vector<Component> Component) : mComponents(Component)
    {
        AddDefaults();
        CheckComponents();
    }

    ///@}
    ///@name Operations
    ///@{

    const PointType& LowerBound() const {
        return Get<PointType>("lower_bound");
    }

    const PointType& UpperBound() const {
        return Get<PointType>("upper_bound");
    }

    const Vector3i& Order() const {
        return Get<Vector3i>("polynomial_order");
    }

    const IntegrationMethodType& IntegrationMethod() const {
        return Get<IntegrationMethodType>("integration_method");
    }

    const Vector3i& NumberOfElements() const {
        return Get<Vector3i>("number_of_elements");;
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

    /// @brief Set values to given component. If parameter is not stored yet, parameter is added.
    /// @tparam type
    /// @param rName Name (Key) of parameter.
    /// @param rValue New value of parameter.
    template<typename type>
    void Set(const std::string& rName, const type& rValue){
        auto p_value = pFind<type>(rName);

        if( p_value ){ // Component already exists -> Override.
            *p_value = rValue;
        }
        else { // Component does not exist -> Add new one.
            mComponents.push_back(Component(rName, rValue));
        }

        CheckComponents();
    }

    /// @brief Get parameter value. Throws an error if parameter is not available.
    /// @tparam type
    /// @param rName Name (Key) of parameter.
    /// @return const type&
    template<typename type>
    const type& Get( const std::string& rName ) const {
        const auto p_value = pFind<type>(rName);
        if( p_value ){
            return *p_value;
        }

        std::string error_message = "Parameter :: Get :: Component: '" + rName + "' not found.\n";
        throw std::runtime_error(error_message);
    }

    /// @brief Print all paramters (Name, Value).
    /// @param rOStream
    void PrintInfo(std::ostream& rOStream) const {
        rOStream << "Parameters: \n";
        for( auto& value : mComponents ){
            rOStream << value.Name() << ": ";
            std::visit(PrintVisit(rOStream), value.Get());
            rOStream << "\n";
        }
    }

private:

    ///@}
    ///@name Private Operations
    ///@{

    /// @brief Adds all defaults values (mDefaultComponents) to current set of parameters.
    void AddDefaults(){
        for( auto& value : mDefaultComponents){
            AddValueIfTypesMatch<PointType>(value);
            AddValueIfTypesMatch<Vector3i>(value);
            AddValueIfTypesMatch<bool>(value);
            AddValueIfTypesMatch<double>(value);
            AddValueIfTypesMatch<IndexType>(value);
            AddValueIfTypesMatch<std::string>(value);
            AddValueIfTypesMatch<IntegrationMethodType>(value);
        }
    }

    /// @brief Adds Value to mComponents if Types of 'type' and 'rValues' match.
    /// @tparam type
    /// @param rValues New Component.
    template<typename type>
    void AddValueIfTypesMatch(const Component& rValues){
        auto p_type_id = std::visit(TypeVisit{}, rValues.Get());
        if( *p_type_id == typeid(type) ){
            auto p_value = pFind<type>(rValues.Name());
            if( !p_value ){
                mComponents.push_back( rValues );
            }
        }
    }

    /// @brief Checks if all components are part of mAllAvailableComponents. Also checks if all
    ///        components have the correct Types.
    void CheckComponents() const {
        for( auto& r_components : mComponents ){
            const auto& current_name = r_components.Name();

            const auto p_pair_found = std::find_if( mAllAvailableComponents.begin(), mAllAvailableComponents.end(),
                [&current_name](const auto& rPair) { return rPair.first == current_name; } );

            if( p_pair_found != mAllAvailableComponents.end() ){
                const auto p_current_type_info = std::visit(TypeVisit{}, r_components.Get());
                const auto p_ref_type_info = p_pair_found->second;
                if( !(*p_current_type_info ==  *p_ref_type_info) ){ // If type is wrong.
                    std::string error_message = "Parameters :: CheckComponents :: Name: '" + current_name +
                            "' is not provided with correct Type.\n";
                    throw std::runtime_error(error_message);
                }
            }
            else { // If name is part of mAllAvailableComponents.
                std::string error_message = "Parameters :: CheckComponents :: Name: '" + current_name +
                    "' is not a valid Parameter.\n";
                throw std::runtime_error(error_message);
            }
        }
    }

    /// @brief Returns Value of parameter if it exists. Returns nullptr otherwise. Const version.
    /// @tparam type
    /// @param rName Name (Key) of parameter.
    /// @return const type*
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

    /// @brief Returns Value of parameter if it exists. Returns nullptr otherwise. Non-Const version.
    /// @tparam type
    /// @param rName Name (Key) of parameter.
    /// @return const type*
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

    ///@}
    ///@name Private Member Variables
    ///@{

    std::vector<Component> mComponents{};
    static const std::vector<Component> mDefaultComponents;
    static const std::vector<std::pair<std::string, const std::type_info*>> mAllAvailableComponents;
    ///@}

}; // End class Parameters
///@} End TIBRA classes

std::ostream& operator<< (std::ostream& rOStream, const Parameters& rThis);

} // End namespace tibra

#endif // PARAMETERS_INCLUDE_H
