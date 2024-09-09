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
#include <cxxabi.h>
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
 * @brief  Stores Key/Value pairs. Keys are stored as TVariantKey. Values are stored as VariantType, which is hard-coded to:
 *         std::<PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType>.
 * @tparam TVariantKey. This should be an enum class, wrapped inside an std::variant<enum class>.
**/
template<typename TVariantKey>
class KeyValuePair {

public:
    ///@name Type definitions
    ///@{

    typedef std::variant<PointType, Vector3i, bool, double, unsigned long, std::string, IntegrationMethodType> VariantType;

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor
    /// @param NewKey Key. Is stored as TVariantKey.
    /// @param NewKeyName std:::string (Is stored to be able to print Key)
    /// @param NewValue Value. Is stored as std::variant<VariantType>
    KeyValuePair(TVariantKey NewKey, std::string NewKeyName, VariantType NewValue, bool Set )
        : mKey(NewKey), mKeyName(NewKeyName), mValue(NewValue), mSet(Set)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns Key
    /// @return const std::variant&
    const TVariantKey& GetKey() const {
        return mKey;
    }

    /// @brief Returns Value (no exception version)
    /// @return const std::variant&
    template<typename TValueType>
    const TValueType& GetValueNoCheck() const {
        assert(("KeyValuePair :: Value is not set.", mSet));
        const TValueType* p_value = std::get_if<TValueType>(&mValue);
        assert(("KeyValuePair :: Given Value type does not match stored Value type.", p_value != 0));
        return *p_value;
    }

    /// @brief Returns Value
    /// @return const std::variant&
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

    /// @brief Sets Value
    /// @tparam Type of value.
    /// @param NewValue
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

    /// @brief Returns name of key.
    /// @return const std::string&
    const std::string& GetKeyName() const {
        return mKeyName;
    }

    /// @brief Returns true is Value is set.
    /// @return bool
    bool IsSet() const {
        return mSet;
    }

    /// Returns underlying type name of stored value as char*.
    const char* GetValueTypeName() const {
        return std::visit( GetTypeNameVisit{}, mValue);
    }

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

    /// Returns type name of TType as char*
    template<typename TType>
    const char* GetTypeName() const {
        int status;
        return abi::__cxa_demangle(typeid(TType).name(),0,0,&status);
    }

    /// Visit struct to print Type names.
    struct GetTypeNameVisit {
        const char* operator()(const PointType& rValue){return "PointType/std::array<double, 3>"; };
        const char* operator()(const Vector3i& rValue){return "Vector3i/std::array<std::size_t, 3>"; };
        const char* operator()(const unsigned long& rValue){return "std::size_t"; };
        const char* operator()(const double& rValue){return "double"; };
        const char* operator()(const std::string& rValue){return "std::string"; };
        const char* operator()(const bool& rValue){return "bool"; };
        const char* operator()(const IntegrationMethodType& rValue){return "IntegrationMethod"; };
    };

    struct PrintVisit {
        PrintVisit(std::ostream& rOStream) : mOstream(rOStream){}
        void operator()(const PointType& rValue){mOstream << rValue; };
        void operator()(const Vector3i& rValue){mOstream << rValue;};
        void operator()(const unsigned long& rValue){ mOstream << rValue;};
        void operator()(const double& rValue){mOstream << rValue;};
        void operator()(const std::string& rValue){mOstream << rValue;};
        void operator()(const bool& rValue){mOstream << rValue; };
        void operator()(const IntegrationMethodType& rValue){ mOstream << rValue; };

    private:
        std::ostream& mOstream;
    };

    ///@}
    ///@name Private member variables
    ///@{

    TVariantKey mKey;
    std::string mKeyName;
    VariantType mValue;
    bool mSet;
    ///@}
};

/**
 * @class  Dictionary.
 * @author Manuel Messmer
 * @brief Stores json-like data sets. Each key points either to the next level in the dictionary or to a base value type, e.g. double.
 *        On the lowest level, KeyValuePair's are stored. @see KeyValuePair. Possible ValueTypes are: PointType, Vector3i, bool, double, unsigned long, std::string, IntegrationMethodType.
 *        The KeyTypes are supposed to be enum classes and can be passed as a template parameter pack: ...TEnumKeys.
 * @tparam TVariantKey. This should be an enum class, wrapped inside an std::variant<enum class>.
**/
template<typename... TEnumKeys>
class Dictionary
{
public:

    ///@name Type definitions
    ///@{

    typedef std::variant<TEnumKeys...> TVariantKey;

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor
    /// @tparam TKey of this dictionary level.
    /// @param rKey Key of this dictionary level.
    template<typename TKeyType>
    Dictionary(TKeyType rKey, std::string rKeyName) : mKeyType(rKey), mKeyName(rKeyName) {
        mKeyTypeName = GetTypeName<TKeyType>();
    }

    /// Destructor
    virtual ~Dictionary() = default;
    /// Copy Constructor
    Dictionary(const Dictionary& ) = delete;
    /// Assignement operator
    Dictionary& operator=(const Dictionary& ) = delete;
    /// Move constructor
    Dictionary(Dictionary&& ) = default;
    /// Move assignement operator
    Dictionary& operator=(Dictionary&& ) = default;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns const ref to sub dictionary (no check version).
    /// @tparam TKey.
    /// @param QueryKey
    /// @return const Dictionary&
    template<typename TKeyType>
    const Dictionary& GetSubDictionaryNoCheck(TKeyType QueryKey) const {
        assert("Given Key type does not match stored Key type: " && std::get_if<TKeyType>(&mSubDictionaryDummyKey) != 0  );
        IndexType index = static_cast<IndexType>(QueryKey);
        return *mSubDictionaries[index];
    }

    /// @brief Returns const ref to next dictionery level.
    /// @tparam TKey.
    /// @param QueryKey
    /// @return const Dictionary&
    template<typename TKeyType>
    const Dictionary& operator [](TKeyType QueryKey) const {
        if ( std::get_if<TKeyType>(&mSubDictionaryDummyKey) ){
            IndexType index = static_cast<IndexType>(QueryKey);
            return *mSubDictionaries[index];
        }
        QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mSubDictionaryKeyTypeName << ").\n";
    }

    /// @brief Returns ref to dictionery corresponding to given key.
    /// @tparam TKey.
    /// @param QueryKey
    /// @return Dictionary&
    template<typename TKeyType>
    Dictionary& operator [](TKeyType QueryKey) {
        if ( std::get_if<TKeyType>(&mSubDictionaryDummyKey) ){
            IndexType index = static_cast<IndexType>(QueryKey);
            return *mSubDictionaries[index];
        }
        QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mSubDictionaryKeyTypeName << ").\n";
    }

    /// @brief Returns value to given Key (no exception version).
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param QueryKey
    /// @return const TValueType&
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValueNoCheck(TKeyType QueryKey) const {
        assert( ("Given Key type does not match stored Key type.", std::get_if<TKeyType>(&mDataDummyKey) ) != 0 );
        IndexType index = static_cast<IndexType>(QueryKey);
        return mData[index].template GetValueNoCheck<TValueType>();
    }

    /// @brief Returns value to given Key.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param QueryKey
    /// @return const TValueType&
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValue(TKeyType QueryKey) const {
        if( std::get_if<TKeyType>(&mDataDummyKey) ) {
            IndexType index = static_cast<IndexType>(QueryKey);
            return mData[index].template GetValue<TValueType>();
        }
        QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mDataKeyTypeName << ").\n";
    }

    /// @brief Sets Value to given Key.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param QueryKey
    /// @param NewValue
    template<typename TKeyType, typename TValueType>
    void SetValue(TKeyType QueryKey, TValueType NewValue) {
        if( std::get_if<TKeyType>(&mDataDummyKey) ){
            IndexType index = static_cast<IndexType>(QueryKey);
            mData[index].SetValue(NewValue);
        } else {
            QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mDataKeyTypeName << ").\n";
        }
    }


    template<typename TKeyType>
    bool IsSet(TKeyType QueryKey){
        if( std::get_if<TKeyType>(&mDataDummyKey) ){
            IndexType index = static_cast<IndexType>(QueryKey);
            return mData[index].IsSet();
        } else {
            QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mDataKeyTypeName << ").\n";
        }
    }

    /// @brief Add values to the current level of the dictionary. Should only be called in derived class to define the default dictionary.
    /// @tparam TKey Type of keys. All keys within one level must have the same type (same enum).
    /// @tparam ...TValues Parameter pack for types.
    /// @param NewEntries
    template<typename TKeyType, typename... TValueTypes>
    void AddValues(std::tuple<std::tuple<TKeyType, std::string, TValueTypes, bool>...> NewEntries) {
        mDataDummyKey = std::get<0>(std::get<0>(NewEntries));
        mDataKeyTypeName = GetTypeName<TKeyType>();
        UnpackTuple(NewEntries);
    }

    /// @brief Returns numver of sub-dictionaries stored.
    /// @return IndexType
    const IndexType NumberOfSubDictionaries() const {
        return mSubDictionaries.size();
    }

    /// @brief Adds an empty AddEmptySubDictionary
    /// @tparam TKeyType
    /// @param rKey
    /// @param rKeyName
    /// @return Dictionary&
    template<typename TKeyType>
    Dictionary& AddEmptySubDictionary(TKeyType rKey, std::string rKeyName) {
        if( mSubDictionaries.size() == 0 ) {
            mSubDictionaryDummyKey = rKey;
            mSubDictionaryKeyTypeName = GetTypeName<TKeyType>();
        }
        else {
            for(const auto& level : mSubDictionaries) {
                if( !std::get_if<TKeyType>(&level->GetKey() ) ) {
                    QuESo_ERROR << "The keys on the same level must be of same type. Given Key type: '" << GetTypeName<TKeyType>()
                        << ". Stored Key type on this level: '" << level->GetKeyTypeName() << "'.\n";
                }
            }
        }
        // Add new dictionary
        mSubDictionaries.push_back( MakeUnique<Dictionary>(rKey, rKeyName) );
        return *mSubDictionaries.back();
    }

    /// @brief Prints this dictionary
    /// @param rOStream
    /// @param Indent
    void PrintInfo(std::ostream& rOStream, std::string& Indent) const {
        rOStream << Indent << mKeyName << " {\n";
        std::string new_indent = Indent + "   ";
        for( const auto& key_value_pair : mData ){
            rOStream << new_indent << key_value_pair.GetKeyName() << ": ";
            key_value_pair.PrintValue(rOStream);
            rOStream << '\n';
        }
        for( const Unique<Dictionary<TEnumKeys...>>& p_sub_dictionary : mSubDictionaries  ) {
            p_sub_dictionary->PrintInfo(rOStream, new_indent);
        }
        if( Indent == "" ){
            rOStream << "}";
        } else {
            rOStream << Indent << "}\n";
        }
    }

    friend std::ostream& operator<< (std::ostream& rOStream, const Dictionary<TEnumKeys ...>& rThis) {
        std::string indent = "";
        rThis.PrintInfo(rOStream, indent);
        return rOStream;
    }

private:

    ///@}
    ///@name Struct and Function definitions to extract tuples.
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
        UnpackTuple<I+1>(Tuple);
    }

    ///@}
    ///@name Private operations.
    ///@{

 /// @brief Returns key of current dictionary level
    /// @return TVariantKey
    const TVariantKey& GetKey() const {
        return mKeyType;
    }

    /// @brief Returns name of Key type.
    /// @return const std::string&
    const std::string& GetKeyTypeName() const {
        return mKeyTypeName;
    }

    /// Returns type name of TType as char*
    template<typename TType>
    const char* GetTypeName() const {
        int status;
        return abi::__cxa_demangle(typeid(TType).name(),0,0,&status);
    }

    ///@}
    ///@name Private member variables.
    ///@{

    //// Key related members. (Key of this dictionary)
    TVariantKey mKeyType;
    std::string mKeyName;
    std::string mKeyTypeName;

    //// Current Value related members. (Values stored in this dictionary)
    TVariantKey mDataDummyKey;
    std::string mDataKeyTypeName;
    std::vector<KeyValuePair<TVariantKey>> mData;

    //// Members related to sub-dictionaries.
    // Dummy key to check if a given key (e.g. given to a functions) has the correct key type.
    // Since the key is used as index for mSubDictionaries, it has to be checked before accesing mSubDictionaries.
    TVariantKey mSubDictionaryDummyKey;
    std::string mSubDictionaryKeyTypeName;
    std::vector<Unique<Dictionary<TEnumKeys...>>> mSubDictionaries;

    ///@}
};

///@}
} // End queso namespace.

#endif // End CONDITION_INCLUDE_HPP