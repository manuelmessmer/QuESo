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

#ifndef DICTIONARY_INCLUDE_HPP
#define DICTIONARY_INCLUDE_HPP

//// STL includes
#include <iostream>
#include <variant>
#include <cassert>
#include <set>
#include <vector>
#include <type_traits>

/// Project includes
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  KeyValuePair
 * @author Manuel Messmer
 * @brief  Stores Key/Value pairs. Keys are stored as TVariantKeyType. Values are stored as VariantValueType, which is hard-coded to:
 *         std::variant<PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType, GridTypeType>.
 * @tparam TVariantKey. This should be an pack of enum classes, wrapped inside an std::variant<...enum class>.
 * @see    Dictionary, which uses KeyValuePair.
**/
template<typename TVariantKeyType>
class KeyValuePair {

public:
    ///@name Type definitions
    ///@{

    typedef std::variant<PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType, GridTypeType> VariantValueType;

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor.
    /// @param NewKey Key. Is stored as TVariantKey.
    /// @param NewKeyName std:::string (Is stored to be able to print Key).
    /// @param NewValue Value. Is stored as VariantValueType/
    /// @param Set True, if value should be set.
    KeyValuePair(TVariantKeyType NewKey, std::string NewKeyName, VariantValueType NewValue, bool Set )
        : mKey(NewKey), mKeyName(NewKeyName), mValue(NewValue), mSet(Set)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns Key
    /// @return const TVariantKeyType&
    const TVariantKeyType& GetKey() const {
        return mKey;
    }

    /// @brief Returns Value (no exception version). Very fast, but does not check types.
    /// @return const TValueType&.
    /// @see GetValue(). <- This version is slower, as it checks the given types.
    template<typename TValueType>
    const TValueType& GetValueNoCheck() const {
        assert(("KeyValuePair :: Value is not set.", mSet));
        const TValueType* p_value = std::get_if<TValueType>(&mValue);
        assert(("KeyValuePair :: Given Value type does not match stored Value type.", p_value != 0));
        return *p_value;
    }

    /// @brief Returns Value (exception version). Slightly slower as GetValueNoCheck(), as it checks the given types.
    ///        On C++ level GetValueNoCheck() should be used.
    /// @tparam TValueType.
    /// @return const TValueType&.
    /// @see GetValue(). <- This version is faster, as it does not check the given types.
    template<typename TValueType>
    const TValueType& GetValue() const {
        if ( mSet ){
            const TValueType* p_type = std::get_if<TValueType>(&mValue);
            if( p_type ) {
                return *p_type;
            }
            QuESo_ERROR << "For key: '" << GetKeyName() << "' - Given Value type (" << GetTypeName<TValueType>() << ") does not match stored Value type ("
                << GetValueTypeName() << ").\n";
        } else {
           QuESo_ERROR << "The Value corresponding to the given Key (" << mKeyName << ") is not set.\n";
        }
    }

    /// @brief Sets Value.
    /// @tparam TValueType of value.
    /// @param NewValue.
    /// @see SetValueWithAmbiguousType() <- This version cast ambiguous types, if possible, e.g. allows to cast 0 -> 0.0.
    template<typename TValueType>
    void SetValue(TValueType NewValue) {
        if( std::get_if<TValueType>(&mValue) ){
            mSet = true;
            mValue = NewValue;
        } else {
            QuESo_ERROR << "For key: '" << GetKeyName() << "' - Given Value type (" << GetTypeName<TValueType>() << ") does not match stored Value type ("
                << GetValueTypeName() << ").\n";
        }
    }

    /// @brief Sets Value and cast 0 -> 0.0 and [0, 0, 0] -> [0.0, 0.0, 0.0] if possible.
    /// @tparam TValueType of value.
    /// @param NewValue
    /// @see SetValue().
    template<typename TValueType>
    void SetValueWithAmbiguousType(TValueType NewValue) {
        bool is_set = false;
        if constexpr ( std::is_same<TValueType, IndexType>::value ) { // If incoming value is IndexType
            if( std::get_if<double>(&mValue) ){                       // If stored value is double
                // Cast 0 -> 0.0
                SetValue( static_cast<double>(NewValue) );
                is_set = true;
            }
        }
        else if constexpr ( std::is_same<TValueType, Vector3i>::value ) { // If incoming value is Vector3i
            if( std::get_if<Vector3d>(&mValue) ){                         // If stored value is Vector3d
                // Cast [0, 0, 0] -> [0.0, 0.0, 0.0]
                SetValue( Vector3d{ static_cast<double>(NewValue[0]),
                                    static_cast<double>(NewValue[1]),
                                    static_cast<double>(NewValue[2])} );
                is_set = true;
            }
        }
        if( !is_set ){
            SetValue(NewValue);
        }
    }

    /// @brief Returns name of key.
    /// @return const std::string&
    const std::string& GetKeyName() const {
        return mKeyName;
    }

    /// @brief Returns true if Value is set.
    /// @return bool
    bool IsSet() const {
        return mSet;
    }

    /// @brief Returns underlying type name of stored value as const char*.
    /// @return const char*
    const char* GetValueTypeName() const {
        return std::visit( GetTypeNameVisit{}, mValue);
    }

    /// @brief Prints Value
    /// @param rOstream
    void PrintValue(std::ostream& rOstream) const {
        if( mSet ){
            std::visit( PrintVisit(rOstream), mValue );
        } else {
            rOstream << "Not Set.";
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
        const char* operator()(const PointType& rValue){return "PointType/std::array<double, 3>"; };
        const char* operator()(const Vector3i& rValue){return "Vector3i/std::array<std::size_t, 3>"; };
        const char* operator()(const IndexType& rValue){return "std::size_t"; };
        const char* operator()(const double& rValue){return "double"; };
        const char* operator()(const std::string& rValue){return "std::string"; };
        const char* operator()(const bool& rValue){return "bool"; };
        const char* operator()(const IntegrationMethodType& rValue){return "IntegrationMethod"; };
        const char* operator()(const GridTypeType& rValue){return "GridType"; };
    };

    /// Visit struct to print values.
    struct PrintVisit {

        PrintVisit(std::ostream& rOStream) : mOstream(rOStream){}
        void operator()(const PointType& rValue){mOstream << rValue; };
        void operator()(const Vector3i& rValue){mOstream << rValue;};
        void operator()(const IndexType& rValue){ mOstream << rValue;};
        void operator()(const double& rValue){mOstream << rValue;};
        void operator()(const std::string& rValue){mOstream << rValue;};
        void operator()(const bool& rValue){std::string out = (rValue) ? "true" : "false"; mOstream << out; };
        void operator()(const IntegrationMethodType& rValue){ mOstream << rValue; };
        void operator()(const GridTypeType& rValue){ mOstream << rValue; };

    private:
        std::ostream& mOstream;
    };

    ///@}
    ///@name Private member variables
    ///@{

    TVariantKeyType mKey;
    std::string mKeyName;
    VariantValueType mValue;
    bool mSet;
    ///@}
};

/**
 * @class  Dictionary.
 * @author Manuel Messmer
 * @brief Stores json-like data sets. Each key points either to a subdictionary or to a base value type, e.g. double.
 *        Key/Value-paires are stored as 'KeyValuePair'. The KeyTypes are supposed to be enum classes and can be passed as a template
 *        parameter pack: ...TEnumKeys.
 *        The main idea is to derive from Dictionary and to set the respective Keys and default Values.
 * @see KeyValuePair. Possible ValueTypes are: PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType, GridTypeType.
 * @see Settings. Settings derives from Dictionary.
 * @tparam TEnumKeys. This should be a pack of enum classes.
**/
template<typename... TEnumKeys>
class Dictionary
{
public:

    ///@name Type definitions
    ///@{

    typedef std::variant<TEnumKeys...> TVariantKey; // Each stpred key can take these values.

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor
    /// @tparam TKeyType of this dictionary level.
    /// @param rKey Key of this dictionary level.
    /// @param rKeyName Name of the Key. Used to print the Key name.
    template<typename TKeyType>
    Dictionary(TKeyType rKey, std::string rKeyName) :
        mKey(rKey), mKeyName(rKeyName), mKeyTypeName(GetTypeName<TKeyType>())
    {
    }

    /// Destructor
    virtual ~Dictionary() = default;
    /// Copy Constructor
    Dictionary(const Dictionary& rDict) = default;
    /// Assignement operator
    Dictionary& operator=(const Dictionary& rDict) = delete;
    /// Move constructor
    Dictionary(Dictionary&& rDict) = default;
    /// Move assignement operator
    Dictionary& operator=(Dictionary&& rDict) = default;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns const ref to sub dictionary (no exception version).
    ///        This version is very fast, as it does not check the given type. This should be used on C++ level.
    /// @tparam TKeyType.
    /// @param QueryKey
    /// @return const Dictionary&
    /// @see operator [] <- Version with exception handling.
    template<typename TKeyType>
    const Dictionary& GetSubDictionaryNoCheck(TKeyType QueryKey) const {
        assert("Given Key type does not match stored Key type: " && std::get_if<TKeyType>(&mSubDictionaryDummyKey) != 0  );
        const IndexType index = static_cast<IndexType>(QueryKey);
        return mSubDictionaries[index];
    }

    /// @brief Returns ref to subdictionary corresponding to given key (exception version, const version).
    ///        This version is slower than GetSubDictionaryNoCheck(), as it checks the given type.
    /// @tparam TKeyType.
    /// @param QueryKey
    /// @return const Dictionary&
    /// @see GetSubDictionaryNoCheck() <- Fast version without exception handling. Should be used, if access Keys are hard-coded.
    template<typename TKeyType>
    const Dictionary& operator [](TKeyType QueryKey) const {
        if constexpr(std::is_same_v<TKeyType, int>
                  || std::is_same_v<TKeyType, unsigned int>
                  || std::is_same_v<TKeyType, unsigned long>) { // Accept different integer types.
            if ( std::get_if<IndexType>(&mSubDictionaryDummyKey) ){
                QuESo_ERROR_IF( QueryKey < 0 ) << "Key must be non-negative.\n";
                const IndexType index = static_cast<IndexType>(QueryKey);
                return mSubDictionaries[index];
            }
        } else {
            if ( std::get_if<TKeyType>(&mSubDictionaryDummyKey) ){
                const IndexType index = static_cast<IndexType>(QueryKey);
                return mSubDictionaries[index];
            }
        }
        QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mSubDictionaryKeyTypeName << ").\n";
    }

    /// @brief Returns ref to subdictionary corresponding to given key (non-const version).
    /// @tparam TKeyType.
    /// @param QueryKey
    /// @return Dictionary&
    /// @see GetSubDictionaryNoCheck() <- Fast version without exception handling. Should be used, if access Keys are hard-coded.
    template<typename TKeyType>
    Dictionary& operator [](TKeyType QueryKey) {
        if constexpr(std::is_same_v<TKeyType, int>
                  || std::is_same_v<TKeyType, unsigned int>
                  || std::is_same_v<TKeyType, unsigned long>) { // Accept different integer types.
            if ( std::get_if<IndexType>(&mSubDictionaryDummyKey) ){
                QuESo_ERROR_IF( QueryKey < 0 ) << "Key must be non-negative.\n";
                const IndexType index = static_cast<IndexType>(QueryKey);
                return mSubDictionaries[index];
            }
        } else {
            if ( std::get_if<TKeyType>(&mSubDictionaryDummyKey) ){
                const IndexType index = static_cast<IndexType>(QueryKey);
                return mSubDictionaries[index];
            }
        }
        QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mSubDictionaryKeyTypeName << ").\n";
    }

    /// @brief Returns value to given Key (no exception version).
    ///        This version is very fast, as it does not check the given type. This should be used on C++ level.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param QueryKey
    /// @return const TValueType&
    /// @see GetValue() <- Version with exception handling.
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValueNoCheck(TKeyType QueryKey) const {
        assert( ("Given Key type does not match stored Key type.", std::get_if<TKeyType>(&mDataDummyKey) ) != 0 );
        const IndexType index = static_cast<IndexType>(QueryKey);
        return mData[index].template GetValueNoCheck<TValueType>();
    }

    /// @brief Returns value to given Key (exception version).
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param QueryKey
    /// @return const TValueType&
    /// @see GetValueNoCheck() <- Version without exception handling. Should be used, if access Keys are hard-coded.
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValue(TKeyType QueryKey) const {
        if( std::get_if<TKeyType>(&mDataDummyKey) ) {
            const IndexType index = static_cast<IndexType>(QueryKey);
            return mData[index].template GetValue<TValueType>();
        }
        QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mDataKeyTypeName << ").\n";
    }

    /// @brief Sets Value to given Key.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param QueryKey
    /// @param NewValue
    /// @see SetValueWithAmbiguousType() <- allows to cast ambiguous types, e.g., casts 0 -> 0.0, if possible.
    template<typename TKeyType, typename TValueType>
    void SetValue(TKeyType QueryKey, TValueType NewValue) {
        if( std::get_if<TKeyType>(&mDataDummyKey) ){
            const IndexType index = static_cast<IndexType>(QueryKey);
            if constexpr(std::is_same_v<TValueType, int>
                      || std::is_same_v<TValueType, unsigned int>
                      || std::is_same_v<TValueType, unsigned long>) { // Accept different integer types.
                QuESo_ERROR_IF(NewValue < 0) << "Value must be non-negative.\n";
                mData[index].SetValue(static_cast<IndexType>(NewValue));
            } else {
                mData[index].SetValue(NewValue);
            }
        } else {
            QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mDataKeyTypeName << ").\n";
        }
    }

    /// @brief Sets Value to given Key. Casts 0 -> 0.0 and [0, 0, 0] -> [0.0, 0.0, 0.0] if possible.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param QueryKey
    /// @param NewValue
    /// @param SetValue() <- Version without type casting.
    template<typename TKeyType, typename TValueType>
    void SetValueWithAmbiguousType(TKeyType QueryKey, TValueType NewValue) {
        if( std::get_if<TKeyType>(&mDataDummyKey) ){
            const IndexType index = static_cast<IndexType>(QueryKey);
            if constexpr(std::is_same_v<TValueType, int>
                      || std::is_same_v<TValueType, unsigned int>
                      || std::is_same_v<TValueType, unsigned long>) { // Accept different integer types.
                QuESo_ERROR_IF(NewValue < 0) << "Value must be non-negative.\n";
                mData[index].SetValueWithAmbiguousType(static_cast<IndexType>(NewValue));
            } else {
                mData[index].SetValueWithAmbiguousType(NewValue);
            }
        } else {
            QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mDataKeyTypeName << ").\n";
        }
    }

    /// @brief Returns true if Value to given Key is set.
    /// @tparam TKeyType.
    /// @param QueryKey
    /// @return bool
    template<typename TKeyType>
    bool IsSet(TKeyType QueryKey) const {
        if( std::get_if<TKeyType>(&mDataDummyKey) ){
            const IndexType index = static_cast<IndexType>(QueryKey);
            return mData[index].IsSet();
        }
        QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mDataKeyTypeName << ").\n";
    }

    /// @brief Adds values to the current level of the dictionary. Should only be called in derived class to define the default dictionary.
    ///        Layout: std::tuple( std::tuple(Enum Key, string KeyName, ValueType Value, bool IsSet ) )
    ///        If IsSet=false, 'Value' is only used to deduce the associated type.
    /// @tparam TKeyType of keys. All keys within one level must have the same type (same enum).
    /// @tparam ...TValueTypes Parameter pack for types.
    /// @param NewEntries
    template<typename TKeyType, typename... TValueTypes>
    void AddValues(std::tuple<std::tuple<TKeyType, std::string, TValueTypes, bool>...> NewEntries) {
        mDataDummyKey = std::get<0>(std::get<0>(NewEntries)); // We store one key, to be able to test for the Key Type.
        mDataKeyTypeName = GetTypeName<TKeyType>();
        UnpackTuple(NewEntries);
    }

    /// @brief Returns number of subdictionaries stored.
    /// @return IndexType
    IndexType NumberOfSubDictionaries() const {
        return mSubDictionaries.size();
    }

    /// @brief Adds an empty subdictionary.
    /// @tparam TKeyType
    /// @param rKey
    /// @param rKeyName
    /// @return Dictionary&
    template<typename TKeyType>
    Dictionary& AddEmptySubDictionary(TKeyType rKey, std::string rKeyName) {
        if( mSubDictionaries.size() == 0 ) {
            mSubDictionaryDummyKey = rKey; // Store key to be able to test the type later.
            mSubDictionaryKeyTypeName = GetTypeName<TKeyType>();
        }
        else {
            for(const auto& sub_dict : mSubDictionaries) {
                if( !std::get_if<TKeyType>(&sub_dict.GetKey() ) ) {
                    QuESo_ERROR << "The keys on the same level must be of same type. Given Key type: '" << GetTypeName<TKeyType>()
                        << ". Stored Key type on this level: '" << sub_dict.GetKeyTypeName() << "'.\n";
                }
            }
        }
        // Add new dictionary
        mSubDictionaries.push_back( Dictionary(rKey, rKeyName) );
        return mSubDictionaries.back();
    }

    /// @brief Throws error if any stored value on this level is not set.
    void CheckIfValuesAreSet() const {
        for( const auto& r_key_val_pair : mData  ){
            QuESo_ERROR_IF( !r_key_val_pair.IsSet() ) << "Variable: " << r_key_val_pair.GetKeyName() << " is not set.\n";
        }
    }

    /// @brief Prints this dictionary
    /// @param rOStream
    /// @param Indent Should be "" when called. This variable is used for recursive function calls.
    void PrintInfo(std::ostream& rOStream, std::string& Indent) const {
        rOStream << Indent << mKeyName << " {\n";
        std::string new_indent = Indent + "   ";
        for( const auto& key_value_pair : mData ){
            rOStream << new_indent << key_value_pair.GetKeyName() << ": ";
            key_value_pair.PrintValue(rOStream);
            rOStream << '\n';
        }
        for( const Dictionary<TEnumKeys...>& r_sub_dictionary : mSubDictionaries  ) {
            r_sub_dictionary.PrintInfo(rOStream, new_indent);
        }
        if( Indent == "" ){
            rOStream << "}";
        } else {
            rOStream << Indent << "}\n";
        }
    }

    // Return begin iterator to list of subdictionaries
    typename std::vector<Dictionary<TEnumKeys...>>::iterator begin_sub_dicts() {
        return mSubDictionaries.begin();
    }

    // Return end iterator to list of subdictionaries
    typename std::vector<Dictionary<TEnumKeys...>>::iterator end_sub_dicts() {
        return mSubDictionaries.end();
    }

    // Operator<<
    friend std::ostream& operator<< (std::ostream& rOStream, const Dictionary<TEnumKeys ...>& rThis) {
        std::string indent = "";
        rThis.PrintInfo(rOStream, indent);
        return rOStream;
    }

private:

    ///@}
    ///@name Struct and function definitions to extract tuples.
    ///@{

    template <typename>
    struct is_tuple: std::false_type {};

    template <typename ...T>
    struct is_tuple<std::tuple<T...>>: std::true_type {};

    // Function to iterate through all values of the tuple. I equals number of values in tuple.
    template <size_t I = 0, typename TKeyType, typename... TValueTypes>
    typename std::enable_if<I == sizeof...(TValueTypes), void>::type
    UnpackTuple(std::tuple<std::tuple<TKeyType, std::string, TValueTypes, bool>...> Tuple) {
        return; // If iterated through all values of tuple, then simply return.
    }

    // Function to iterate through all values of the tuple. I equals number of values in tuple.
    template <size_t I = 0, typename TKeyType, typename... TValueTypes>
    typename std::enable_if<(I < sizeof...(TValueTypes)), void>::type
    UnpackTuple(std::tuple<std::tuple<TKeyType, std::string, TValueTypes, bool>...> Tuple) {
        auto& r_tuple = std::get<I>(Tuple);
        if(  static_cast<IndexType>(std::get<0>(r_tuple)) != I ){
            QuESo_ERROR << "Tuples are not given in correct order. \n";
        }

        mData.push_back(KeyValuePair(TVariantKey(std::get<0>(r_tuple)), std::get<1>(r_tuple), std::get<2>(r_tuple), std::get<3>(r_tuple)));
        UnpackTuple<I+1>(Tuple); // Recursive call.
    }

    ///@}
    ///@name Private operations.
    ///@{

    /// @brief Returns Key of current dictionary level
    /// @return  const TVariantKey&
    const TVariantKey& GetKey() const {
        return mKey;
    }

    /// @brief Returns name of Key type.
    /// @return const std::string&
    const std::string& GetKeyTypeName() const {
        return mKeyTypeName;
    }

    /// Returns type name of TType as char*
    template<typename TType>
    const char* GetTypeName() const {
        // Names are not really pretty, but abi::__cxa__demangle only works for gcc.
        return typeid(TType).name();
    }

    ///@}
    ///@name Private member variables.
    ///@{

    //// Key related members. (Key of this dictionary)
    TVariantKey mKey;
    std::string mKeyName;
    std::string mKeyTypeName;

    //// Current Value related members. (Values stored in this dictionary).
    // This key (mDataDummyKey) allows to test if any given Key type matches the stored Key types.
    // It is usually the Key of the first entry of mData.
    TVariantKey mDataDummyKey;
    std::string mDataKeyTypeName;
    std::vector<KeyValuePair<TVariantKey>> mData;

    //// Members related to subdictionaries.
    // This key (mSubDictionaryDummyKey) allows to test if any given Key type matches the stored Key types.
    // Since the key is used as index for mSubDictionaries, it has to be checked before accesing mSubDictionaries.
    TVariantKey mSubDictionaryDummyKey;
    std::string mSubDictionaryKeyTypeName;
    std::vector<Dictionary<TEnumKeys...>> mSubDictionaries;

    ///@}
};

///@} // End QuESo Classes
} // End queso namespace.

#endif // End DICTIONARY_INCLUDE_HPP