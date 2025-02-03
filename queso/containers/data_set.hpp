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

#ifndef DATA_SET_INCLUDE_HPP
#define DATA_SET_INCLUDE_HPP

//// STL includes
#include <variant>

/// Project includes
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  DataSet
 * @author Manuel Messmer
 * @brief  Stores Key/Value pairs. Keys are stored as TVariantKeyType. Values are stored as VariantValueType, which is hard-coded to:
 *         std::variant<PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType, GridTypeType>.
 * @tparam TVariantKey. This should be an pack of enum classes, wrapped inside an std::variant<...enum class>.
 * @see    Dictionary, which uses DataSet.
**/
class DataSet {
public:
    ///@name Type definitions
    ///@{

    typedef std::variant<std::monostate, PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType, GridTypeType> VariantValueType;

    // Helper tag type
    template<class TKeyType>
    struct TypeTag {
        typedef TKeyType KeyType;
    };

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor.
    /// @param NewKey Key. Is stored as TVariantKeyType.
    /// @param NewKeyName std::string (Is stored to be able to print Key).
    /// @param NewValue Value. Is stored as VariantValueType.
    /// @param Set True, if value should be set.
    /// @note Even if DataSet is not set, NewValue must be given such that the underlying type
    ///       (which this DataSet is supposed to hold) can be deduced.
    ///       When later mValue is actually set, its type is checked against the type of the already stored dummy value.
    ///       This scenario can not be handled with std::monostate.
    template<typename TTypeTag>
    DataSet(TTypeTag) {
        typedef decltype(key::GetKeyBaseType<typename TTypeTag::KeyType>()) BaseType;
        static_assert( std::is_same<typename BaseType::KeyToWhat, key::KeyToDataSet>::value );
        mpKeyInformation = MakeUnique<typename BaseType::KeyInfo>();
        mData.resize(mpKeyInformation->GetNumberOfKeys());
    }

    /// @brief Sets Value to given Key provided as Enum.
    /// @tparam TKeyType
    /// @tparam TValueType
    /// @param QueryKey (Enum)
    /// @param rNewValue
    /// @see Only throws in DEBUG mode.
    template<typename TKeyType, typename TValueType>
    void SetValue(TKeyType QueryKey, const TValueType& rNewValue,
                  std::enable_if_t<std::is_enum<TKeyType>::value>* = nullptr ) noexcept(NOTDEBUG) {
        
        QuESo_ASSERT( key::IsCorrectType(mpKeyInformation, QueryKey), "Wrong Key Type.\n" );
        QuESo_ASSERT( variable::IsCorrectType<TValueType>(QueryKey), "Wrong Value Type.\n" );
        
        const IndexType index = static_cast<IndexType>(QueryKey);
        if constexpr(std::is_same_v<TValueType, int>
                  || std::is_same_v<TValueType, unsigned int>
                  || std::is_same_v<TValueType, unsigned long>) { // Accept different integer types.
            if constexpr (std::is_same_v<TValueType, int>) {
                QuESo_ASSERT(rNewValue >= 0, "Value must be non-negative.\n");
            }
            mData[index] = static_cast<IndexType>(rNewValue);
        } else {
            mData[index] = rNewValue;
        }
    }

    /// @brief Sets Value to given Key (given as string).
    /// @tparam TValueType
    /// @param rQueryKeyName
    /// @param rNewValue
    /// @see SetValueWithAmbiguousType() <- allows to cast ambiguous types, e.g., casts 0 -> 0.0, if possible.
    /// @note Should only be used in Python. There is also a version that takes an actual KeyType instead of a string.
    template<typename TValueType>
    void SetValue(const std::string& rQueryKeyName, const TValueType& rNewValue) {
        IndexType index = mpKeyInformation->GetKeyValue(rQueryKeyName);
        
        if constexpr(std::is_same_v<TValueType, int>
                || std::is_same_v<TValueType, unsigned int>
                || std::is_same_v<TValueType, unsigned long>) { // Accept different integer types.
            if constexpr (std::is_same_v<TValueType, int>) {
                QuESo_ERROR_IF(rNewValue < 0) << "Value must be non-negative.\n";
            }
             mData[index] = static_cast<IndexType>(rNewValue);
        } else {
            mData[index] = rNewValue;
        }
    }

    /// @brief Returns the Value to the given Key provided as Enum.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param QueryKey (Enum)
    /// @return const TValueType&
    /// @note Only throws in DEBUG mode.
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValue(TKeyType QueryKey,    
                               std::enable_if_t<std::is_enum<TKeyType>::value>* = nullptr ) const noexcept(NOTDEBUG){
        QuESo_ASSERT( key::IsCorrectType(mpKeyInformation, QueryKey), "Given Key Type does not match stored Key type.\n" );
        const int index = static_cast<int>(QueryKey);
        QuESo_ASSERT( std::get_if<std::monostate>(&mData[index]) == 0, "Value is not set.");
        const TValueType* p_value = std::get_if<TValueType>(&mData[index]);
        QuESo_ASSERT(p_value != 0, "Given Value type does not match stored Value type.");
        return *p_value;
    }

private:

    ///@}
    ///@name Private operations
    ///@{

    /// Returns type name of TType as const char*
    template<typename TType>
    const char* GetTypeName() const {
        // Names are not really pretty, but abi::__cxa__demangle only works for gcc.
        return typeid(TType).name();
    }

    /// Visit struct to get Type names.
    struct GetTypeNameVisit {
        std::string operator()(const PointType& rValue){return "PointType/std::array<double, 3>"; };
        std::string operator()(const Vector3i& rValue){return "Vector3i/std::array<std::size_t, 3>"; };
        std::string operator()(const IndexType& rValue){return "std::size_t"; };
        std::string operator()(const double& rValue){return "double"; };
        std::string operator()(const std::string& rValue){return "std::string"; };
        std::string operator()(const bool& rValue){return "bool"; };
        std::string operator()(const IntegrationMethodType& rValue){return "IntegrationMethod"; };
        std::string operator()(const GridTypeType& rValue){return "GridType"; };
    };

    /// Visit struct to print values in JSON format.
    struct PrintVisit {

        PrintVisit(std::ostream& rOStream) : mOstream(rOStream){}
        void operator()(const PointType& rValue){mOstream << '[' << rValue[0] << ", " << rValue[1] << ", " << rValue[2] << ']'; };
        void operator()(const Vector3i& rValue){mOstream << '[' << rValue[0] << ", " << rValue[1] << ", " << rValue[2] << ']';};
        void operator()(const IndexType& rValue){ mOstream << rValue;};
        void operator()(const double& rValue){mOstream << rValue;};
        void operator()(const std::string& rValue){mOstream << '\"' << rValue << '\"';};
        void operator()(const bool& rValue){std::string out = (rValue) ? "true" : "false"; mOstream << out; };
        void operator()(const IntegrationMethodType& rValue){ mOstream << '\"' << rValue << '\"'; };
        void operator()(const GridTypeType& rValue){ mOstream << '\"' << rValue << '\"'; };

    private:
        std::ostream& mOstream;
    };

    ///@}
    ///@name Private struct and function definitions to extract tuples.
    ///@{



    ///@}
    ///@name Private member variables
    ///@{

    Unique<key::KeyInformation> mpKeyInformation;
    std::vector<VariantValueType> mData;
    ///@}
}; // End class DataSet

} // End queso namespace

#endif // End DATA_SET_INCLUDE_HPP