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
 * @brief  Stores values that can be accessed with Keys from a KeySet. The size of the DataSet depends on the size of the KeySet.
 *         The KeySet must be registered in "includes/register_keys.hpp".
 *         The respective values are stored in a std::vector<VariantValueType> instance. Possible VariantValueTypes are:
 *         PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType, GridTypeType.
**/
class DataSet {
public:
    ///@name Type definitions
    ///@{

    typedef std::variant<std::monostate, PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType, GridTypeType> VariantValueType;

    /// Helper to check for TType is a scoped int enum.
    template<typename TType>
    using is_scoped_int_enum = std::integral_constant<bool, !std::is_convertible<TType,int>::value &&
                                                             std::is_enum<TType>::value>;


    // Helper to check if TType can be cast to std::size_t.
    template<typename TType, typename = void>
    struct is_castable_to_size_t : std::false_type {};

    template<typename TType>
    struct is_castable_to_size_t<TType, std::void_t<decltype(static_cast<std::size_t>(std::declval<TType>()))>> : std::true_type {};


    // Helper to check if TType is_index (allows signed integers, but no bool and char types).
    template<typename TType>
    using is_index = std::integral_constant<bool, is_castable_to_size_t<TType>::value &&
                                                  std::is_integral<TType>::value &&
                                                  !std::is_same<TType, bool>::value &&
                                                  !std::is_same<TType, char>::value &&
                                                  !std::is_same<TType, unsigned char>::value &&
                                                  !std::is_same<TType, signed char>::value >;

    // Helper to check if TType is_index (allows signed integers, but no bool and char types).
    template<typename TType>
    static constexpr bool is_index_v = is_index<TType>::value;

    // Helper to check if TType is_unsigned_index (only accepts unsigned integers, also no bool and char types).
    template<typename TType>
    using is_unsigned_index = std::integral_constant<bool, is_index_v<TType> &&
                                                           std::is_unsigned<TType>::value>;

    // Helper to check if TType is_unsigned_index (only accepts unsigned integers, also no bool and char types).
    template<typename TType>
    static constexpr bool is_unsigned_index_v = is_unsigned_index<TType>::value;

    ///@}
    ///@name Life cycle
    ///@{

    // Helper type tag to templetize the constructor of DataSet.
    template<class TKeySetInfoType>
    struct TypeTag {
        typedef TKeySetInfoType KeySetInfoType;
    };

    /// @brief Constructor
    /// @tparam TTypeTag
    /// @param TTypeTag contains typedef to KeySetInfoType: TTypeTag<KeySetInfoType>.
    template<typename TTypeTag>
    DataSet(TTypeTag) : mpKeySetInfo( MakeUnique<typename TTypeTag::KeySetInfoType>() ) {
        static_assert( std::is_same<typename TTypeTag::KeySetInfoType::KeyToWhat, key::KeyToValue>::value );
        mData.resize(mpKeySetInfo->GetNumberOfKeys());
    }

    /// @brief Sets the value to the given Key.
    /// @tparam TKeyType
    /// @tparam TValueType
    /// @param rQueryKey
    /// @param rNewValue
    /// @note Only throws in DEBUG mode.
    template<typename TKeyType, typename TValueType>
    void SetValue(const TKeyType& rQueryKey, const TValueType& rNewValue,
                  std::enable_if_t<is_scoped_int_enum<typename TKeyType::KeyValueType>::value>* = nullptr ) noexcept(NOTDEBUG) {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeyToWhat, key::KeyToValue>::value );
        static_assert( std::is_same<typename TKeyType::KeyToWhat, TValueType>::value ||
                       (std::is_same<typename TKeyType::KeyToWhat, IndexType>::value &&
                        is_unsigned_index<TValueType>::value) );

        QuESo_ASSERT( mpKeySetInfo->IsCorrectKeyType(rQueryKey), "Given Key: '" + rQueryKey.Name() + "' is of wrong type.\n" );
        if constexpr( is_unsigned_index_v<TValueType> ) { // Accept different unsigned integer types.
            // is_unsigned_index_v implies that rNewValue is castable to IndexType.
            mData[rQueryKey.Index()] = static_cast<IndexType>(rNewValue);
        } else {
            mData[rQueryKey.Index()] = rNewValue;
        }
    }

    /// @brief Sets the value to the given Key (given as string).
    /// @tparam TValueType
    /// @param rQueryKeyName
    /// @param rNewValue
    /// @note Should only be used in Python. There is also a version that takes an actual KeyType instead of a string.
    template<typename TValueType>
    void SetValue(const std::string& rQueryKeyName, const TValueType& rNewValue) {
        const auto p_key = mpKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName
                               << "'. Possible names are: " + mpKeySetInfo->GetAllKeyNames() << "\n.";

        const auto target_type_index = p_key->VariableTypeIndex();

        // Case 1: Handling TValueTypes that are convertible to 'IndexType'. (signed and unsigned types are allowed.)
        if constexpr( is_index_v<TValueType> ) {
            if( target_type_index == std::type_index(typeid(double)) ){ // Target type: 'double'.
                // Allow cast from IndexType to double.
                mData[p_key->Index()] = static_cast<double>(rNewValue);
                return;
            }
            else if( target_type_index == std::type_index(typeid(IndexType)) ) { // Target type: 'IndexType'.
                if constexpr( std::is_signed_v<TValueType> ) { // Now, TValueType must be unsigned.
                    QuESo_ERROR_IF(rNewValue < 0) << "Value for Key: '" << p_key->Name() << "' must be non-negative.\n";
                }
                mData[p_key->Index()] = static_cast<IndexType>(rNewValue);
                return;
            }
            QuESo_ERROR << "The given key: '" << p_key->Name() << "' accesses a variable of type: '" << p_key->VariableTypeName()
                        << "'. However, the given value type is: 'std::size_t / int'.\n";
        }
        // Case 2: TValueType is Vector3i.
        else if constexpr( std::is_same<TValueType, Vector3i>::value ) {
            if( target_type_index == std::type_index(typeid(Vector3d)) ){ // Target type is 'Vector3d'.
                // Allow cast from Vector3i to Vector3d.
                mData[p_key->Index()] = Vector3d{ static_cast<double>(rNewValue[0]),
                                                  static_cast<double>(rNewValue[1]),
                                                  static_cast<double>(rNewValue[2]) };
                return;
            }
            else if ( target_type_index == std::type_index(typeid(Vector3i)) ) { // Target type is 'Vector3i'.
                mData[p_key->Index()] = rNewValue;
                return;
            }
            QuESo_ERROR << "The given key: '" << p_key->Name() << "' accesses a variable of type: '" << p_key->VariableTypeName()
                        << "'. However, the given value type is: 'Vector3i / std::array<double,3>'.\n";
        }
        else {
            // Case 3: Target type and given type match.
            if( target_type_index == std::type_index(typeid(TValueType)) ){
                mData[p_key->Index()] = rNewValue;
                return;
            }
            QuESo_ERROR << "The given key: '" << p_key->Name() << "' accesses a variable of type: '" << p_key->VariableTypeName()
                        << "'. However, the given value type is: '" << key::KeyToValue::template Visit<TValueType>::visit() << "'.\n";
        }
    }



    /// @brief Returns the value to the given Key (fast version, does not throw).
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return const TValueType&
    /// @note Only throws in debug mode. In realese it does not check if the value is set and if TKeyType is of correct type.
    /// @see GetValue(): Thorws if value is not set.
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValueFast(const TKeyType& rQueryKey,
                                   std::enable_if_t<is_scoped_int_enum<typename TKeyType::KeyValueType>::value>* = nullptr ) const noexcept(NOTDEBUG) {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeyToWhat, key::KeyToValue>::value );
        static_assert( std::is_same<typename TKeyType::KeyToWhat, TValueType>::value );

        QuESo_ASSERT( mpKeySetInfo->IsCorrectKeyType(rQueryKey), "Given Key: '" + rQueryKey.Name() + "' is of wrong type.\n" );
        QuESo_ASSERT( !std::get_if<std::monostate>(&mData[rQueryKey.Index()]),  "Value to Key: '" + rQueryKey.Name() + "' is not set.\n");

        return std::get<TValueType>(mData[rQueryKey.Index()]);
    }

    /// @brief Returns the value to the given Key.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return const TValueType&
    /// @note In release mode it throws if the given key is not set. In debug mode it also throws if TKeyType is of wrong type.
    ///       TKeyType should be known and therefore always be correct.
    /// @see GetValueFast(): Never throws in release.
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValue(const TKeyType& rQueryKey,
                               std::enable_if_t<is_scoped_int_enum<typename TKeyType::KeyValueType>::value>* = nullptr ) const {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeyToWhat, key::KeyToValue>::value );
        static_assert( std::is_same<typename TKeyType::KeyToWhat, TValueType>::value );

        QuESo_ASSERT( mpKeySetInfo->IsCorrectKeyType(rQueryKey), "Given Key: '" + rQueryKey.Name() + "' is of wrong type.\n" );
        QuESo_ERROR_IF( std::get_if<std::monostate>(&mData[rQueryKey.Index()]) ) << "Value to Key: '" + rQueryKey.Name() + "' is not set.\n";

        return std::get<TValueType>(mData[rQueryKey.Index()]);
    }


    /// @brief Returns the value the given Key (given as std::string).
    /// @tparam TValueType
    /// @param rQueryKeyName
    /// @return const TValueType&
    /// @note Should only be used in Python. There is also a version that takes an actual KeyType instead of a string.
    template<typename TValueType>
    const TValueType& GetValue(const std::string& rQueryKeyName) const {
        const auto p_key = mpKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName << "'. Possible names are: " + mpKeySetInfo->GetAllKeyNames();

        QuESo_ERROR_IF( std::get_if<std::monostate>(&mData[p_key->Index()]) ) << "Value to Key: '" + p_key->Name() + "' is not set.\n";
        QuESo_ERROR_IF( p_key->VariableTypeIndex() != std::type_index(typeid(TValueType)) ) << "The given key: '" << p_key->Name()
            << "' accesses a variable of type: '" << p_key->VariableTypeName()
            << "'. However, the given value type is: '" << key::KeyToValue::template Visit<TValueType>::visit() << "'.\n";

        return std::get<TValueType>(mData[p_key->Index()]);
    }

    /// @brief Returns true if the value to the given Key is set.
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return bool
    /// @note Only throws in debug mode.
    template<typename TKeyType>
    bool IsSet(const TKeyType& rQueryKey,
               std::enable_if_t<is_scoped_int_enum<typename TKeyType::KeyValueType>::value>* = nullptr ) const noexcept(NOTDEBUG) {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeyToWhat, key::KeyToValue>::value );
        QuESo_ASSERT( mpKeySetInfo->IsCorrectKeyType(rQueryKey), "Given Key: '" + rQueryKey.Name() + "' is of wrong type.\n" );
        return !std::get_if<std::monostate>(&mData[rQueryKey.Index()]);
    }

    /// @brief Returns true if the value to the given Key (given as std::string) is set.
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return bool
    /// @note Should only be used in Python. There is also a version that takes an actual KeyType instead of a string.
    bool IsSet(const std::string& rQueryKeyName) const {
        const auto p_key = mpKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName << "'. Possible names are: " + mpKeySetInfo->GetAllKeyNames();

        return !std::get_if<std::monostate>(&mData[p_key->Index()]);
    }

    /// @brief Returns true if given key (given as std::string) exists.
    /// @param rQueryKeyName
    /// @return bool
    /// @note Should only be used in Python. In C++ the correct key types must be known.
    bool Has(const std::string& rQueryKeyName) const {
        return mpKeySetInfo->pGetKey(rQueryKeyName) != nullptr;
    }

private:

    ///@}
    ///@name Private member variables
    ///@{

    Unique<key::KeySetInfo> mpKeySetInfo;
    std::vector<VariantValueType> mData;

    ///@}
}; // End class DataSet

} // End queso namespace

#endif // End DATA_SET_INCLUDE_HPP