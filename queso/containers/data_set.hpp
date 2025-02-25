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
    template<class TKeySetInfoType>
    struct TypeTag {
        typedef TKeySetInfoType KeySetInfoType;
    };

    /// Helper to check for scoped int enums
    template<typename TType>
    using is_scoped_int_enum = std::integral_constant<bool, !std::is_convertible<TType,int>::value &&
                                                             std::is_enum<TType>::value>;


    // Helper to check if TType can be cast to std::size_t
    template<typename TType, typename = void>
    struct is_castable_to_size_t : std::false_type {};

    template<typename TType>
    struct is_castable_to_size_t<TType, std::void_t<decltype(static_cast<std::size_t>(std::declval<TType>()))>> : std::true_type {};


    // Helper to check if TType is_index (allows signed integers, but no bool and char types)
    template<typename TType>
    using is_index = std::integral_constant<bool, is_castable_to_size_t<TType>::value &&
                                                  std::is_integral<TType>::value &&
                                                  !std::is_same<TType, bool>::value &&
                                                  !std::is_same<TType, char>::value &&
                                                  !std::is_same<TType, unsigned char>::value &&
                                                  !std::is_same<TType, signed char>::value >;

    // Helper to check if TType is_index (allows signed integers, but no bool and char types)
    template<typename TType>
    static constexpr bool is_index_v = is_index<TType>::value;

    // Helper to check if TType is_index (only accepts unsigned integers, also no bool and char types)
    template<typename TType>
    using is_unsigned_index = std::integral_constant<bool, is_index_v<TType> &&
                                                           std::is_unsigned<TType>::value>;

    // Helper to check if TType is_index (only accepts unsigned integers, also no bool and char types)
    template<typename TType>
    static constexpr bool is_unsigned_index_v = is_unsigned_index<TType>::value;

    ///@}
    ///@name Life cycle
    ///@{

    template<typename TTypeTag>
    DataSet(TTypeTag) {
        static_assert( std::is_same<typename TTypeTag::KeySetInfoType::KeyToWhat, key::KeyToValue>::value );
        mpKeySetInfo = MakeUnique<typename TTypeTag::KeySetInfoType>();
        mData.resize(mpKeySetInfo->GetNumberOfKeys());
    }

    /// @brief Sets Value to given Key.
    /// @tparam TKeyType
    /// @tparam TValueType
    /// @param QueryKey (Enum)
    /// @param rNewValue
    /// @see Only throws in DEBUG mode.
    template<typename TKeyType, typename TValueType>
    void SetValue(const TKeyType& rQueryKey, const TValueType& rNewValue,
                  std::enable_if_t<is_scoped_int_enum<typename TKeyType::KeyValueType>::value>* = nullptr ) noexcept(NOTDEBUG) {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeyToWhat, key::KeyToValue>::value );
        static_assert( std::is_same<typename TKeyType::KeyToWhat, TValueType>::value ||
                       (std::is_same<typename TKeyType::KeyToWhat, IndexType>::value &&
                        is_unsigned_index<TValueType>::value) );

        QuESo_ASSERT( mpKeySetInfo->IsCorrectKeyType(rQueryKey), "Given Key is of wrong type.\n" );
        if constexpr( is_unsigned_index_v<TValueType> ) { // Accept different integer types.
            mData[rQueryKey.Index()] = static_cast<IndexType>(rNewValue);
        } else {
            mData[rQueryKey.Index()] = rNewValue;
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
        // This will throw if the key does not exist.
        const auto p_key = mpKeySetInfo->pGetKey(rQueryKeyName);

        if constexpr( is_index_v<TValueType> ) { // Given type : convertable to 'IndexType' (Can be signed).
            if( p_key->VariableTypeIndex() == std::type_index(typeid(double)) ){ // Target type: 'double'.
                // Allow cast from double to IndexType.
                mData[p_key->Index()] = static_cast<double>(rNewValue);
                return;
            }
            else if( p_key->VariableTypeIndex() == std::type_index(typeid(IndexType)) ) { // Target type: 'IndexType'.
                if constexpr (std::is_signed<TValueType>::value) {
                    QuESo_ERROR_IF(rNewValue < 0) << "Value for Key: '" << p_key->Name() << "' must be non-negative.\n";
                }
                mData[p_key->Index()] = static_cast<IndexType>(rNewValue);
                return;
            }
            QuESo_ERROR << "The given key accesses: '" << p_key->VariableTypeName() << "'. However, the given value is: 'std::size_t or int'.\n";
        } else if constexpr( std::is_same<TValueType, Vector3i>::value ) { // Given type is Vector3i.
            if( p_key->VariableTypeIndex() == std::type_index(typeid(Vector3d)) ){ // Target type is 'Vector3d'.
                // Allow cast from Vector3i to Vector3d.
                mData[p_key->Index()] = Vector3d{ static_cast<double>(rNewValue[0]),
                                                  static_cast<double>(rNewValue[1]),
                                                  static_cast<double>(rNewValue[2]) };
                return;
            }
            else if ( p_key->VariableTypeIndex() == std::type_index(typeid(Vector3i)) ) { // Target type is 'Vector3i'.
                mData[p_key->Index()] = rNewValue;
                return;
            }
            QuESo_ERROR << "The given key accesses: '" << p_key->VariableTypeName() << "'. However, the given value is: 'Vector3i'.\n";
        } else {
            if( p_key->VariableTypeIndex() == std::type_index(typeid(TValueType)) ){ // Given type and target type match.
                mData[p_key->Index()] = rNewValue;
                return;
            }
            QuESo_ERROR << "The given key accesses: '" << p_key->VariableTypeName() << "'. However, the given value is: '"
                << key::KeyToValue::template Visit<TValueType>::visit() << "'.\n";
        }
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

    Unique<key::KeySetInfo> mpKeySetInfo;
    std::vector<VariantValueType> mData;
    ///@}
}; // End class DataSet

} // End queso namespace

#endif // End DATA_SET_INCLUDE_HPP