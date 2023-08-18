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

namespace queso {

///@name  QuESo Globals
///@{

enum IntegrationMethod {Gauss, Gauss_Reduced1, Gauss_Reduced2, GGQ_Optimal, GGQ_Reduced1, GGQ_Reduced2};
typedef enum IntegrationMethod IntegrationMethodType;

enum ConditionType {Neumann, Dirichlet};
typedef enum ConditionType ConditionTypeType;

///@}
///@name  QuESo Classes
///@{

/**
 * @class  ParamCondition (Base class)
 * @author Manuel Messmer
 * @brief  Interface for ParamConditionNeumann and ParamConditionDirichlet.
**/
class ParamCondition {
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ParamCondition(IndexType Id, const std::string& rFilename, const Vector3d& rPrescribed )
        : mId(Id), mFilename( rFilename ), mPrescribed( rPrescribed )
    {
    }
    /// Destructor
    virtual ~ParamCondition()= default;

    ///@}
    ///@name Operations
    ///@{

    /// Returns type of condition.
    virtual ConditionTypeType Type() const = 0;

    /// Returns filename.
    virtual const std::string& GetFilename() const {
        return mFilename;
    }

    /// Returns condition id.
    virtual IndexType GetId() const {
        return mId;
    }

    /// Returns prescribed values.
    virtual const Vector3d& GetPrescribed() const {
        return mPrescribed;
    }

    /// Returns penalty factor.
    virtual double GetPenaltyFactor() const {
        QuESo_ERROR("ParamCondition::GetPenaltyFactor") << "Calling base class. Penalty factor is only available for 'ParamConditionDirichlet'.\n";
    }

private:
    ///@}
    ///@name Private Members
    ///@{
    IndexType mId;
    std::string mFilename;
    Vector3d mPrescribed;
    ///@}
}; // End class ParamCondition

/**
 * @class  ParamConditionNeumann
 * @author Manuel Messmer
 * @brief  Container for neumann condition related parameters.
**/
class ParamConditionNeumann : public ParamCondition {
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ParamConditionNeumann(IndexType Id, const std::string& rFilename, const Vector3d& rPrescribed)
        : ParamCondition(Id, rFilename, rPrescribed )
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /// Returns type of condition.
    ConditionTypeType Type() const override {
        return ConditionType::Neumann;
    }
    ///@}

}; /// End of class ParamConditionNeumann

/**
 * @class  ParamConditionDirichlet
 * @author Manuel Messmer
 * @brief  Container for dirichlet condition related parameters.
**/
class ParamConditionDirichlet : public ParamCondition {
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ParamConditionDirichlet(IndexType Id, const std::string& rFilename, const Vector3d& rPrescribed, double PenaltyFactor)
        : ParamCondition(Id, rFilename, rPrescribed ), mPenaltyFactor(PenaltyFactor)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /// Returns type of condition.
    ConditionTypeType Type() const override {
        return ConditionType::Dirichlet;
    }

    /// Returns penalty factor.
    double GetPenaltyFactor() const override {
        return mPenaltyFactor;
    }

private:
    ///@}
    ///@name Private members
    ///@{
    double mPenaltyFactor;
    ///@}
}; // End of ParamConditionDirichlet class.


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
    typedef std::variant<PointType, Vector3i, bool, double, unsigned long, std::string, IntegrationMethodType> ComponentType;

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
    ///@name Typedefs
    ///@{
    typedef std::vector<Shared<ParamCondition>> ConditionPtrVectorType;
    ///@}

    ///@name std::visit Structs
    ///@{
    struct TypeVisit {
        const std::type_info* operator()(const PointType& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const Vector3i& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const unsigned long& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const double& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const std::string& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const bool& rValue){return &typeid(rValue); };
        const std::type_info* operator()(const IntegrationMethodType& rValue){return &typeid(rValue); };
    };

    struct PrintVisit {
        PrintVisit(std::ostream& rOStream) : mOStream(rOStream){}
        void operator()(const PointType& rValue){mOStream << rValue; };
        void operator()(const Vector3i& rValue){mOStream << rValue;};
        void operator()(const unsigned long& rValue){ mOStream << rValue;};
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
        EnsureValidValues();
        CheckComponents();
    }

    ///@}
    ///@name Operations
    ///@{

    const PointType& LowerBoundXYZ() const {
        return Get<PointType>("lower_bound_xyz");
    }

    const PointType& UpperBoundXYZ() const {
        return Get<PointType>("upper_bound_xyz");
    }

    const PointType& LowerBoundUVW() const {
        return Get<PointType>("lower_bound_uvw");
    }

    const PointType& UpperBoundUVW() const {
        return Get<PointType>("upper_bound_uvw");
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
        return static_cast<IndexType>( Get<unsigned long>("min_num_boundary_triangles") );
    }

    double MomentFittingResidual() const {
        return Get<double>("moment_fitting_residual");
    }

    IndexType EchoLevel() const {
        return static_cast<IndexType>( Get<unsigned long>("echo_level") );
    }

    bool UseCustomizedTrimmedPositions() const{
        /// Only used for testing.
        return Get<bool>("use_customized_trimmed_points");
    }

    IndexType GetPointDistributionFactor() const {
        return static_cast<IndexType>( Get<unsigned long>("init_point_distribution_factor") );
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

        EnsureValidValues();
        CheckComponents();
    }

    /// @brief Returns true if parameter component exists.
    /// @tparam type
    /// @param rName Name (Key) of parameter.
    /// @return bool
    template<typename type>
    bool Has( const std::string& rName ) const {
        const auto p_value = pFind<type>(rName);
        if( p_value ){
            return true;
        }
        return false;
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

        QuESo_ERROR("Parameters::Get") << "Component: '" + rName + "' not found.\n";
    }

    /// @brief Adds neumann condition to parameters.
    /// @param Id Id of condition (Must be unique).
    /// @param rFilename of STL.
    /// @param rPrescribed force.
    void AddNeumannCondition(IndexType Id, const std::string& rFilename, const PointType& rPrescribed ) {
        mConditions.push_back( MakeShared<ParamConditionNeumann>(Id, rFilename, rPrescribed) );
    }

    /// @brief Adds dirichlet condition to parameters.
    /// @param Id Id of condition (Must be unique).
    /// @param rFilename of STL.
    /// @param rPrescribed displacement.
    /// @param PenaltyFactor
    void AddDirichletCondition(IndexType Id, const std::string& rFilename, const PointType& rPrescribed, double PenaltyFactor ){
        mConditions.push_back( MakeShared<ParamConditionDirichlet>(Id, rFilename, rPrescribed, PenaltyFactor) );
    }

    /// @brief Returns number of conditions. Including Neumann and Dirichlet.
    /// @return IndexType
    IndexType NumberOfConditions(){
        return mConditions.size();
    }

    /// @brief Returns filename of condition.
    /// @param Id of condition.
    /// @return const std::string&
    const std::string& GetFilenameOfCondition(IndexType Id) {
        QuESo_ERROR_IF( "Parameters::GetFilenameOfCondition", mConditions[Id]->GetId() != Id ) << "Id does not match index.\n";
        return mConditions[Id]->GetFilename();
    }

    /// @brief Returns vector of ptr to conditions.
    /// @return const ConditionPtrVectorType&
    const ConditionPtrVectorType& GetConditions() const {
        return mConditions;
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
            AddValueIfTypesMatch<unsigned long>(value);
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

    /// @brief Adjust values that are not feasible.
    void EnsureValidValues() {
        // Make sure this value is not numerically zero.
        const double value = Get<double>("min_element_volume_ratio");
        if( value < 0.9e-10 ){
            Set<double>("min_element_volume_ratio", 1e-10);
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
                    QuESo_ERROR("Parameters::CheckComponents") << "Name: '" + current_name +
                        "' is not provided with correct Type.\n";
                }
            }
            else { // If name is part of mAllAvailableComponents.
                QuESo_ERROR("Parameters::CheckComponents") << "Name: '" + current_name +
                    "' is not a valid Parameter.\n";
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
    std::vector<Shared<ParamCondition>> mConditions{};

    inline static const std::vector<Component> mDefaultComponents = {
        Component("echo_level", 0UL),
        Component("embedding_flag", true),
        Component("initial_triangle_edge_length", 1.0),
        Component("min_num_boundary_triangles", 500UL),
        Component("moment_fitting_residual", 1.0e-10),
        Component("min_element_volume_ratio", 1.0e-3),
        Component("b_spline_mesh", true),
        Component("init_point_distribution_factor", 1UL),
        Component("polynomial_order", Vector3i(2UL, 2UL, 2UL) ),
        Component("integration_method", IntegrationMethod::Gauss),
        Component("use_customized_trimmed_points", false) };

    inline static const std::vector<std::pair<std::string, const std::type_info*>> mAllAvailableComponents = {
        std::make_pair<std::string, const std::type_info*>("input_filename", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("postprocess_filename", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("echo_level", &typeid(unsigned long) ),
        std::make_pair<std::string, const std::type_info*>("embedding_flag", &typeid(bool) ),
        std::make_pair<std::string, const std::type_info*>("lower_bound_xyz", &typeid(PointType) ),
        std::make_pair<std::string, const std::type_info*>("upper_bound_xyz", &typeid(PointType) ),
        std::make_pair<std::string, const std::type_info*>("lower_bound_uvw", &typeid(PointType) ),
        std::make_pair<std::string, const std::type_info*>("upper_bound_uvw", &typeid(PointType) ),
        std::make_pair<std::string, const std::type_info*>("b_spline_mesh", &typeid(bool) ),
        std::make_pair<std::string, const std::type_info*>("polynomial_order", &typeid(Vector3i) ),
        std::make_pair<std::string, const std::type_info*>("number_of_elements", &typeid(Vector3i) ),
        std::make_pair<std::string, const std::type_info*>("initial_triangle_edge_length", &typeid(double) ),
        std::make_pair<std::string, const std::type_info*>("min_num_boundary_triangles", &typeid(unsigned long) ),
        std::make_pair<std::string, const std::type_info*>("moment_fitting_residual", &typeid(double) ),
        std::make_pair<std::string, const std::type_info*>("min_element_volume_ratio", &typeid(double) ),
        std::make_pair<std::string, const std::type_info*>("init_point_distribution_factor", &typeid(unsigned long) ),
        std::make_pair<std::string, const std::type_info*>("integration_method", &typeid(IntegrationMethodType) ),
        std::make_pair<std::string, const std::type_info*>("use_customized_trimmed_points", &typeid(bool) ) };

    ///@}

}; // End class Parameters
///@} End QuESo classes

std::ostream& operator<< (std::ostream& rOStream, const Parameters& rThis);

} // End namespace queso

#endif // PARAMETERS_INCLUDE_H
