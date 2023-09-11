// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef PARAMETERS_INCLUDE_H
#define PARAMETERS_INCLUDE_H

////Project includes
#include "includes/define.hpp"
#include "containers/variant_data_container.hpp"

namespace queso {


///@}
///@name  QuESo Classes
///@{

/**
 * @class  Container to store all global parameters.
 * @author Manuel Messmer
**/
class GlobalParameters : public VariantDataContainer {
public:
    typedef VariantDataContainer BaseType;
    using BaseType::Get;
    using BaseType::Set;

    /// Default constructor
    GlobalParameters() : VariantDataContainer() {
        AddDefaults();
        CheckComponents();
    }

    /// Constructor
    GlobalParameters(BaseType::ComponentVectorType Component) : VariantDataContainer(Component)
    {
        AddDefaults();
        EnsureValidValues();
        CheckComponents();
    }

    /// @brief Provides default parameters to base class.
    /// @return const BaseType::ComponentVectorType&
    const BaseType::ComponentVectorType& GetDefaults() const override {
        return mDefaultComponents;
    }

    /// @brief Provides aivalable parameteters to base class.
    /// @return const BaseType::AvailableComponentVectorType&
    const BaseType::AvailableComponentVectorType& GetAvailableComponents() const override {
        return mAllAvailableComponents;
    }

    /// @brief Adjust values that are not feasible.
    void EnsureValidValues() override {
        // Make sure this value is not numerically zero.
        const double value = Get<double>("min_element_volume_ratio");
        if( value < 0.9e-10 ){
            BaseType::Set<double>("min_element_volume_ratio", 1e-10);
        }
    }

private:
    ///@}
    ///@name Private Member Variables
    ///@{

    inline static const BaseType::ComponentVectorType mDefaultComponents = {
        Component("echo_level", 0UL),
        Component("embedding_flag", true),
        Component("initial_triangle_edge_length", 1.0),
        Component("min_num_boundary_triangles", 500UL),
        Component("moment_fitting_residual", 1.0e-10),
        Component("min_element_volume_ratio", 1.0e-3),
        Component("b_spline_mesh", true),
        Component("lower_bound_uvw", PointType(-1.0, -1.0, -1.0)),
        Component("upper_bound_uvw", PointType(1.0, 1.0, 1.0)),
        Component("init_point_distribution_factor", 1UL),
        Component("polynomial_order", Vector3i(2UL, 2UL, 2UL) ),
        Component("integration_method", IntegrationMethod::Gauss),
        Component("use_customized_trimmed_points", false) };

    inline static const BaseType::AvailableComponentVectorType mAllAvailableComponents = {
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
}; // End class GlobalParameters

/**
 * @class  Container to store all condition parameters.
 * @author Manuel Messmer
**/
class ConditionParameters : public VariantDataContainer {
public:
    typedef VariantDataContainer BaseType;
    using BaseType::Get;
    using BaseType::Set;

    /// Default constructor
    ConditionParameters(std::string type) : VariantDataContainer() {
        Set<std::string>("type", type);
        CreateAvailableComponentList();
        AddDefaults();
        CheckComponents();
    }

    /// Constructor
    ConditionParameters(BaseType::ComponentVectorType Component) : VariantDataContainer(Component)
    {
        CreateAvailableComponentList();
        AddDefaults();
        EnsureValidValues();
        CheckComponents();
    }

    void CreateAvailableComponentList() {
        if( Get<std::string>("type") == "PenaltySupportCondition" ) {
            mAvailableComponents.insert(mAvailableComponents.end(), mAvailableComponentsPenalty.begin(), mAvailableComponentsPenalty.end());
        } else if( Get<std::string>("type") == "LagrangeSupportCondition" ) {
            mAvailableComponents.insert(mAvailableComponents.end(), mAvailableComponentsLagrange.begin(), mAvailableComponentsLagrange.end());
        } else if( Get<std::string>("type") == "SurfaceLoadCondition" ) {
            mAvailableComponents.insert(mAvailableComponents.end(), mAvailableComponentsSurfaceLoad.begin(), mAvailableComponentsSurfaceLoad.end());
        } else if( Get<std::string>("type") == "PressureLoadCondition" ) {
            mAvailableComponents.insert(mAvailableComponents.end(), mAvailableComponentsPressureLoad.begin(), mAvailableComponentsPressureLoad.end());
        } else {
            QuESo_ERROR << "Condition type '" << Get<std::string>("type") << "' is not available. Available types are: "
                << "'PenaltySupportCondition', 'LagrangeSupportCondition', 'SurfaceLoadCondition', 'PressureLoadCondition'\n";
        }
    }

    const BaseType::ComponentVectorType& GetDefaults() const override {
        return mDefaultComponents;
    }

    const BaseType::AvailableComponentVectorType& GetAvailableComponents() const override {
        return mAvailableComponents;
    }

    /// @brief Adjust values that are not feasible.
    void EnsureValidValues() override {}

private:
    ///@}
    ///@name Private Member Variables
    ///@{

    inline static const BaseType::ComponentVectorType mDefaultComponents = { };

    BaseType::AvailableComponentVectorType mAvailableComponents = {
        std::make_pair<std::string, const std::type_info*>("type", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("input_filename", &typeid(std::string) ) };

    inline static const BaseType::AvailableComponentVectorType mAvailableComponentsPenalty = {
        std::make_pair<std::string, const std::type_info*>("penalty_factor", &typeid(double) ),
        std::make_pair<std::string, const std::type_info*>("value", &typeid(PointType) ) };

    inline static const BaseType::AvailableComponentVectorType mAvailableComponentsLagrange = {
        std::make_pair<std::string, const std::type_info*>("value", &typeid(PointType) ) };

    inline static const BaseType::AvailableComponentVectorType mAvailableComponentsSurfaceLoad = {
        std::make_pair<std::string, const std::type_info*>("modulus", &typeid(double) ),
        std::make_pair<std::string, const std::type_info*>("direction", &typeid(PointType) ) };

    inline static const BaseType::AvailableComponentVectorType mAvailableComponentsPressureLoad = {
        std::make_pair<std::string, const std::type_info*>("modulus", &typeid(double) ) };
};


/**
 * @class  Parameters
 * @author Manuel Messmer
 * @brief  Dynamic container for all available parameters. Parameters are split into global parameters and a vector of condition parameters:
 *         one for each condition.
**/
class Parameters {
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    Parameters() {
    }

    /// Constructor
    Parameters(const GlobalParameters& rGlobalParameters) : mGlobalParameters(rGlobalParameters)
    {
    }

    /// Constructor
    Parameters(std::vector<Component> Components) : mGlobalParameters(Components)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void AddGlobalParameters(const GlobalParameters& rGlobalParameters) {
        mGlobalParameters = rGlobalParameters;
    }

    void AddConditionParameters(const ConditionParameters& rConditionParameters) {
        mConditionParameters.push_back(rConditionParameters);
    }

    template<typename type>
    const type& Get( const std::string& rName ) const {
        return mGlobalParameters.Get<type>(rName);
    }

    template<typename type>
    void Set(const std::string& rName, const type& rValue) {
        mGlobalParameters.Set(rName, rValue);
    }

    template<typename type>
    bool Has( const std::string& rName ) const {
        return mGlobalParameters.Has<type>(rName);
    }

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

    void AddCondition( const ConditionParameters rConditionParameter ){
        mConditionParameters.push_back(rConditionParameter);
    }

    /// @brief Returns number of conditions. Including Neumann and Dirichlet.
    /// @return IndexType
    IndexType NumberOfConditions(){
        return mConditionParameters.size();
    }

    /// @brief Returns vector of ptr to conditions.
    /// @return std::vector<ConditionParameters>&
    const std::vector<ConditionParameters> & GetConditions() const {
        return mConditionParameters;
    }

    /// @brief Print all paramters (Name, Value).
    /// @param rOStream
    void PrintInfo(std::ostream& rOStream) const {
        rOStream << "Parameters: \n";
        mGlobalParameters.PrintInfo(rOStream);
    }

private:

    ///@}
    ///@name Private Member Variables
    ///@{

    GlobalParameters mGlobalParameters{};
    std::vector<ConditionParameters> mConditionParameters{};

    ///@}

}; // End class Parameters
///@} End QuESo classes

std::ostream& operator<< (std::ostream& rOStream, const Parameters& rThis);

} // End namespace queso

#endif // PARAMETERS_INCLUDE_H
