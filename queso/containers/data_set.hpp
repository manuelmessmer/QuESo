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

/// Forward declaration
template<typename DataSetType>
class DataSetStringAccess;

///@name QuESo classes
///@{

/// @class  DataSet
/// @author Manuel Messmer
/// @tparam TKeySetValuesTypeTag
/// @brief  Static DataSet for fast data access. 'TKeySetValuesTypeTag' must be registered in 'register_keys.hpp'.
/// @details The DataSet allows to access data quickly by ConstexprKeys -> 'register_keys.hpp', also see 'keys.hpp'.
///          The DataSet also provides methods to access the data via std::strings. Given the KeyName a DynamicKey
///          is deduced from KeySetInfo (which is passed to the constructor). The respective KeySetInfo must also
///          be registered in 'register_keys.hpp'.
///          The main idea is to access the data via ConstexprKeys in C++. The std::string access is meant to be used
///          in Python only. Therefore, the respective functions are private and made available via DataSetStringAccess.
///          This is just a safety guard that these functions are not accidentally used in C++.
/// @see 'register_keys.hpp' and 'keys.hpp'.
template<typename TKeySetValuesTypeTag>
class DataSet {
public:
    ///@name Type definitions
    ///@{

    using VariantValueType = typename TKeySetValuesTypeTag::VariantType;

    /// Helper to check if TType is a key.
    template<typename TType>
    static constexpr bool is_key_v = TKeySetValuesTypeTag::template is_key_v<TType>;

    // Helper to check if TType can be cast to std::size_t.
    template<typename TType, typename = void>
    struct is_castable_to_size_t : std::false_type {};

    // Helper to check if TType can be cast to std::size_t.
    template<typename TType>
    struct is_castable_to_size_t<TType, std::void_t<decltype(static_cast<std::size_t>(std::declval<TType>()))>> : std::true_type {};


    // Helper to check if TType is_index (allows signed integers, but no bool and chars).
    template<typename TType>
    using is_index = std::integral_constant<bool, is_castable_to_size_t<TType>::value &&
                                                  std::is_integral<TType>::value &&
                                                  !std::is_same<TType, bool>::value &&
                                                  !std::is_same<TType, char>::value &&
                                                  !std::is_same<TType, unsigned char>::value &&
                                                  !std::is_same<TType, signed char>::value >;

    // Helper to check if TType is_index (allows signed integers, but no bool and chars).
    template<typename TType>
    static constexpr bool is_index_v = is_index<TType>::value;

    // Helper to check if TType is_unsigned_index (only accepts unsigned integers, also no bool and no chars).
    template<typename TType>
    using is_unsigned_index = std::integral_constant<bool, is_index_v<TType> &&
                                                           std::is_unsigned<TType>::value>;

    // Helper to check if TType is_unsigned_index (only accepts unsigned integers, also no bool and no chars).
    template<typename TType>
    static constexpr bool is_unsigned_index_v = is_unsigned_index<TType>::value;

    ///@}
    ///@name Life cycle
    ///@{

    // Helper type tag to use a constructor template.
    template<class TKeySetInfoType>
    struct KeySetInfoTypeTag {
        using KeySetInfoType = TKeySetInfoType;
    };

    /// @brief Constructor
    /// @tparam TKeySetInfoTypeTag
    /// @param TKeySetInfoTypeTag (must be: TKeySetInfoTypeTag<KeySetInfoType>) contains typedef to KeySetInfoType.
    template<typename TKeySetInfoTypeTag>
    DataSet(TKeySetInfoTypeTag) : mpKeySetInfo( MakeUnique<typename TKeySetInfoTypeTag::KeySetInfoType>() ) {
        using KeySetInfoType = typename TKeySetInfoTypeTag::KeySetInfoType;
        static_assert( std::is_same_v<typename KeySetInfoType::KeySetToWhat, TKeySetValuesTypeTag> );
        mData.resize(mpKeySetInfo->GetNumberOfKeys());
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Sets the value associated with the given key.
    /// @tparam TKeyType
    /// @tparam TValueType
    /// @param rQueryKey
    /// @param rNewValue
    /// @note Only throws in DEBUG mode.
    template<typename TKeyType, typename TValueType>
    void SetValue(const TKeyType& rQueryKey, const TValueType& rNewValue,
                  std::enable_if_t<is_key_v<TKeyType>>* = nullptr ) noexcept(NOTDEBUG) {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeySetToWhat, TKeySetValuesTypeTag>::value );
        static_assert( std::is_same<typename TKeyType::KeyToWhat, TValueType>::value ||
                       (std::is_same<typename TKeyType::KeyToWhat, IndexType>::value && is_unsigned_index<TValueType>::value) );

        QuESo_ASSERT( mpKeySetInfo->IsSameKeySet(rQueryKey.KeySetInfoTypeIndex()),
            "Given Key: '" + std::string(rQueryKey.Name()) + "' is of wrong type.\n" );

        if constexpr( is_unsigned_index_v<TValueType> ) { // Accept different unsigned integer types.
            // Note: is_unsigned_index_v implies that rNewValue is castable to IndexType.
            mData[rQueryKey.Index()] = static_cast<IndexType>(rNewValue);
        } else {
            mData[rQueryKey.Index()] = rNewValue;
        }
    }

    /// @brief Returns the value associated with the given Key (fast version, does not throw).
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return const TValueType&
    /// @note Only throws in debug mode. In realese it does not check if the value is set and if TKeyType is of correct type.
    /// @see GetValue(): Throws if value is not set.
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValueFast(const TKeyType& rQueryKey,
                                   std::enable_if_t<is_key_v<TKeyType>>* = nullptr ) const noexcept(NOTDEBUG) {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeySetToWhat, TKeySetValuesTypeTag>::value );
        static_assert( std::is_same<typename TKeyType::KeyToWhat, TValueType>::value );

        QuESo_ASSERT( mpKeySetInfo->IsSameKeySet(rQueryKey.KeySetInfoTypeIndex()),
            "Given Key: '" + std::string(rQueryKey.Name()) + "' is of wrong type.\n" );
        QuESo_ASSERT( !std::holds_alternative<std::monostate>(mData[rQueryKey.Index()]),
            "Value to Key: '" + std::string(rQueryKey.Name()) + "' is not set.\n");

        return std::get<TValueType>(mData[rQueryKey.Index()]);
    }

    /// @brief Returns the value associated with the given Key.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return const TValueType&
    /// @note In release mode it throws if the given key is not set. In debug mode it addionally throws if TKeyType is of wrong type.
    /// @see GetValueFast(): Never throws in release.
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValue(const TKeyType& rQueryKey,
                               std::enable_if_t<is_key_v<TKeyType>>* = nullptr ) const {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeySetToWhat, TKeySetValuesTypeTag>::value );
        static_assert( std::is_same<typename TKeyType::KeyToWhat, TValueType>::value );

        QuESo_ASSERT( mpKeySetInfo->IsSameKeySet(rQueryKey.KeySetInfoTypeIndex()),
            "Given Key: '" + std::string(rQueryKey.Name()) + "' is of wrong type.\n" );
        QuESo_ERROR_IF( std::holds_alternative<std::monostate>(mData[rQueryKey.Index()]) )
            << "Value to Key: '" << rQueryKey.Name() << "' is not set.\n";

        return std::get<TValueType>(mData[rQueryKey.Index()]);
    }

    /// @brief Returns true if the value associated with the given Key is set.
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return bool
    /// @note Only throws in debug mode.
    template<typename TKeyType>
    bool IsSet(const TKeyType& rQueryKey,
               std::enable_if_t<is_key_v<TKeyType>>* = nullptr ) const noexcept(NOTDEBUG) {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeySetToWhat, TKeySetValuesTypeTag>::value );
        QuESo_ASSERT( mpKeySetInfo->IsSameKeySet(rQueryKey.KeySetInfoTypeIndex()),
            "Given Key: '" + std::string(rQueryKey.Name()) + "' is of wrong type.\n" );

        return !std::holds_alternative<std::monostate>(mData[rQueryKey.Index()]);
    }
    ///@}

private:
    ///@name Private operations
    ///@{

    /// The following operations allow to access the data via KeyNames (std::string).
    /// These functions should only be used in Python. Therefore, they are made private here
    /// and can only be used via the friend class DataSetStringAccess.
    friend class DataSetStringAccess<DataSet>;

    /// @brief Sets the value associated with the given Key (given as string).
    /// @tparam TValueType
    /// @param rQueryKeyName
    /// @param rNewValue
    /// @note Should only be used in Python. There is also a version that takes an actual KeyType instead of a std::string.
    template<typename TValueType>
    void SetValue(const std::string& rQueryKeyName, const TValueType& rNewValue) {
        const auto p_key = mpKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName
                               << "'. Possible names are: " << mpKeySetInfo->GetAllKeyNames() << "\n.";

        const auto target_type_index = p_key->TargetTypeIndex();

        // Case 1: TValueType is convertible to 'IndexType' (signed and unsigned types are allowed).
        if constexpr( is_index_v<TValueType> ) {
            if( target_type_index == std::type_index(typeid(double)) ){ // Target type: 'double'.
                // Allow cast from IndexType to double.
                mData[p_key->Index()] = static_cast<double>(rNewValue);
                return;
            }
            else if( target_type_index == std::type_index(typeid(IndexType)) ) { // Target type: 'IndexType'.
                if constexpr( std::is_signed_v<TValueType> ) { // Now, TValueType must be non-negative.
                    QuESo_ERROR_IF(rNewValue < 0) << "Value for Key: '" << p_key->Name() << "' must be non-negative.\n";
                }
                mData[p_key->Index()] = static_cast<IndexType>(rNewValue);
                return;
            }
            QuESo_ERROR << "The given key: '" << p_key->Name() << "' accesses a variable of type: '" << p_key->TargetTypeName()
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
            QuESo_ERROR << "The given key: '" << p_key->Name() << "' accesses a variable of type: '" << p_key->TargetTypeName()
                        << "'. However, the given value type is: 'Vector3i / std::array<double,3>'.\n";
        }
        else {
            // Case 3: Target type and given type match.
            if( target_type_index == std::type_index(typeid(TValueType)) ){
                mData[p_key->Index()] = rNewValue;
                return;
            }
            QuESo_ERROR << "The given key: '" << p_key->Name() << "' accesses a variable of type: '" << p_key->TargetTypeName()
                        << "'. However, the given value type is: '" << TKeySetValuesTypeTag::template GetValueTypeName<TValueType>() << "'.\n";
        }
    }

    /// @brief Returns the value associated with the given Key (given as std::string).
    /// @tparam TValueType
    /// @param rQueryKeyName
    /// @return const TValueType&
    /// @note Should only be used in Python. There is also a version that takes an actual KeyType instead of a std::string.
    template<typename TValueType>
    const TValueType& GetValue(const std::string& rQueryKeyName) const {
        const auto p_key = mpKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName
                               << "'. Possible names are: " << mpKeySetInfo->GetAllKeyNames() << "\n.";

        QuESo_ERROR_IF( std::holds_alternative<std::monostate>(mData[p_key->Index()]) ) << "Value to Key: '" << p_key->Name() << "' is not set.\n";
        QuESo_ERROR_IF( p_key->TargetTypeIndex() != std::type_index(typeid(TValueType)) ) << "The given key: '" << p_key->Name()
            << "' accesses a variable of type: '" << p_key->TargetTypeName()
            << "'. However, the given value type is: '" << TKeySetValuesTypeTag::template GetValueTypeName<TValueType>() << "'.\n";

        return std::get<TValueType>(mData[p_key->Index()]);
    }

    /// @brief Returns true if the value associated with the given Key (given as std::string) is set.
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return bool
    /// @note Should only be used in Python. There is also a version that takes an actual KeyType instead of a std::string.
    bool IsSet(const std::string& rQueryKeyName) const {
        const auto p_key = mpKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName
                               << "'. Possible names are: " << mpKeySetInfo->GetAllKeyNames() << "\n.";

        return !std::holds_alternative<std::monostate>(mData[p_key->Index()]);
    }

    /// @brief Returns true if the given key (given as std::string) exists.
    /// @param rQueryKeyName
    /// @return bool
    /// @note Should only be used in Python. In C++ the correct key types must be known.
    bool Has(const std::string& rQueryKeyName) const {
        return mpKeySetInfo->pGetKey(rQueryKeyName) != nullptr;
    }

    ///@}
    ///@name Private member variables
    ///@{

    Unique<key::detail::KeySetInfo> mpKeySetInfo;
    std::vector<VariantValueType> mData;

    ///@}
}; // End class DataSet

/// @class DataSetStringAccess (FOR USE IN PYTHON ONLY).
/// @brief Allows to access the data of DataSet via KeyNames (std::strings).
///        The respective members are private within DataSet and made public here.
/// @tparam DataSetType
/// @details DataSetStringAccess is a friend of DataSet
template<typename DataSetType>
class DataSetStringAccess {
public:
    ///@name Static operations
    ///@{

    /// @brief Wrapper for DataSet::SetValue.
    /// @tparam TValueType
    /// @param rDataSet
    /// @param rQueryKeyName (std::string)
    /// @param rNewValue
    template<typename TValueType>
    static void SetValue(DataSetType& rDataSet, const std::string& rQueryKeyName, const TValueType& rNewValue) {
        rDataSet.SetValue(rQueryKeyName, rNewValue);
    }

    /// @brief Wrapper for DataSet::GetValue.
    /// @tparam TValueType
    /// @param rDataSet
    /// @param rQueryKeyName (std::string)
    /// @return const TValueType&
    template<typename TValueType>
    static const TValueType& GetValue(const DataSetType& rDataSet, const std::string& rQueryKeyName) {
        return rDataSet.template GetValue<TValueType>(rQueryKeyName);
    }

    /// @brief Wrapper for DataSet::IsSet.
    /// @param rDataSet
    /// @param rQueryKeyName (std::string)
    /// @return bool
    static bool IsSet(const DataSetType& rDataSet, const std::string& rQueryKeyName) {
        return rDataSet.IsSet(rQueryKeyName);
    }

    /// @brief Wrapper for DataSet::Has.
    /// @param rDataSet
    /// @param rQueryKeyName (std::string)
    /// @return bool
    static bool Has(const DataSetType& rDataSet, const std::string& rQueryKeyName) {
        return rDataSet.Has(rQueryKeyName);
    }

    ///@}
}; // End class DataSetStringAccess<>

///@} QuESo classes

} // End queso namespace

#endif // End DATA_SET_INCLUDE_HPP