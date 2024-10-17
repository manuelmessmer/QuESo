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
#include <vector>
#include <algorithm>
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
    /// @param NewKey Key. Is stored as TVariantKeyType.
    /// @param NewKeyName std::string (Is stored to be able to print Key).
    /// @param NewValue Value. Is stored as VariantValueType.
    /// @param Set True, if value should be set.
    /// @note Even if KeyValuePair is not set, NewValue must be given such that the underlying type
    ///       (which this KeyValuePair is supposed to hold) can be deduced.
    ///       When later mValue is actually set, its type is checked against the type of the already stored dummy value.
    ///       This scenario can not be handled with std::monostate.
    KeyValuePair(TVariantKeyType NewKey, std::string NewKeyName, VariantValueType NewValue, bool Set )
        : mKey(NewKey), mKeyName(NewKeyName), mValue(NewValue), mIsSet(Set)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns Key
    /// @return const TVariantKeyType&
    const TVariantKeyType& GetKey() const noexcept {
        return mKey;
    }

    /// @brief Returns Value.
    /// @return const TValueType&.
    /// @note Only throws in DEBUG mode.
    template<typename TValueType>
    const TValueType& GetValueNoCheck() const noexcept(NOTDEBUG) {
        QuESo_ASSERT(mIsSet, "Value is not set.");
        const TValueType* p_value = std::get_if<TValueType>(&mValue);
        QuESo_ASSERT(p_value != 0, "Given Value type does not match stored Value type.");
        return *p_value;
    }

    /// @brief Returns Value.
    /// @tparam TValueType.
    /// @return const TValueType&.
    /// @note Throws in DEBUG and NONBDEUG mode.
    template<typename TValueType>
    const TValueType& GetValue() const {
        if ( mIsSet ){
            const TValueType* p_type = std::get_if<TValueType>(&mValue);
            if( p_type ) {
                return *p_type;
            }
            QuESo_ERROR << "For key: '" << GetKeyName() << "' - Given Value type (" << GetTypeName<TValueType>()
                << ") does not match stored Value type (" << GetValueTypeName() << ").\n";
        }
        QuESo_ERROR << "The Value corresponding to the given Key (" << mKeyName << ") is not set.\n";
    }

    /// @brief Sets Value.
    /// @tparam TValueType of value.
    /// @param NewValue.
    /// @note Only throws in DEBUG mode.
    template<typename TValueType>
    void SetValueNoCheck(TValueType NewValue) noexcept(NOTDEBUG) {
        const TValueType* p_value = std::get_if<TValueType>(&mValue);
        QuESo_ASSERT(p_value != 0, "Given Value type does not match stored Value type.");
        if( p_value ) {
            mIsSet = true;
            mValue = NewValue;
        }
    }

    /// @brief Sets Value.
    /// @tparam TValueType of value.
    /// @param NewValue.
    /// @note Throws in DEBUG and NONBDEUG mode.
    template<typename TValueType>
    void SetValue(TValueType NewValue) {
        if( std::get_if<TValueType>(&mValue) ){
            mIsSet = true;
            mValue = NewValue;
        } else {
            QuESo_ERROR << "For key: '" << GetKeyName() << "' - Given Value type (" << GetTypeName<TValueType>() << ") does not match stored Value type ("
                << GetValueTypeName() << ").\n";
        }
    }

    /// @brief Sets Value and casts 0 -> 0.0 and [0, 0, 0] -> [0.0, 0.0, 0.0] if possible.
    /// @tparam TValueType of value.
    /// @param NewValue
    /// @note This is useful for Python, where the values come from a JSON file.
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
        if( !is_set ){ // Normal set
            SetValue(NewValue);
        }
    }

    /// @brief Returns name of key.
    /// @return const std::string&
    const std::string& GetKeyName() const noexcept {
        return mKeyName;
    }

    /// @brief Returns true if Value is set.
    /// @return bool
    bool IsSet() const noexcept {
        return mIsSet;
    }

    /// @brief Returns underlying type name of stored value as std::string.
    /// @return std::string
    std::string GetValueTypeName() const {
        return std::visit( GetTypeNameVisit{}, mValue);
    }

    /// @brief Prints Value in JSON format
    /// @param rOstream
    void PrintValue(std::ostream& rOstream) const {
        if( mIsSet ){
            std::visit( PrintVisit(rOstream), mValue );
        } else {
            rOstream << "\"Not Set.\"";
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
    ///@name Private member variables
    ///@{

    TVariantKeyType mKey;
    std::string mKeyName;
    VariantValueType mValue;
    bool mIsSet;

    ///@}
};

/// @see Dictionary
struct DictStarts {
    // Must be int. We need to compare to negative values.
    enum Starts : int {start_subdicts = 0, start_values = 0x0100, start_lists=0x0200};
};

/**
 * @class  Dictionary.
 * @author Manuel Messmer
 * @brief Base class to store json-like data sets. We can store KeyValuePair's, Subdictionaries, or Lists.
 *        The Dictionary is optimized to allow extremely fast access to the underlying data. Therefore, the keys are hardcoded.
 *        The KeyTypes are supposed to be enum classes and can be passed as a template parameter pack: ...TEnumKeys.
 *        The main idea is to derive from Dictionary and to 'hardcode' the respective keys/enums.
 *        Thereby, the keys on the same level, must be of the same type (same enum). ( @see e.g. Settings ).
 * @details KeyValuePair, Subdictionaries and Lists are stored in std::vector<>'s. The enum keys are used as the respective indices,
 *          allowing extremely fast access. Consequently, the enum pointing to the first items of KeyValuePair, Subdictionaries, or
 *          Lists must all be equal to 0. However, we want them to be the same EnumType, since they are all on the same level in the
 *          Dictionary. To solve this problem, DictStarts is introduced to shift the respective enum_keys.
 *          For example, an enum key Type may look like:
 *    enum class ExampleKeys {
 *         first_subdict=DictStarts::start_subdicts, second_subdict,
 *         first_key_value_pair=DictStarts::start_values, second_key_value_pair,
 *         first_list=DictStarts::start_lists, second_list }
 *
 *    According to DictStarts, this equals: enum class ExampleKeys {0, 1, 256, 257, 512, 513}
 *    When we now want to e.g. access the first KeyValuePair, we just have to substract 'start_values = 256' from the enum.
 *    Hence, we get {0, 1, 0, 1, 0, 1}.
 *
 * @see KeyValuePair. Possible ValueTypes are: PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType, GridTypeType.
 * @see Settings. Settings derives from Dictionary.
 * @see ModelInfo. ModelInfo derives from Dictionary.
 * @tparam TEnumKeys. This should be a pack of enum classes.
**/
template<typename... TEnumKeys>
class Dictionary
{
public:

    ///@name Type definitions
    ///@{

    /// std::monostate allows to not set the key.
    typedef std::variant<std::monostate, TEnumKeys...> TVariantKey; // Each stored key can take these values.
    typedef KeyValuePair<TVariantKey> KeyValuePairType;

    /// Type trait to check for scoped int enums
    template<typename TType>
    using is_scoped_int_enum = std::integral_constant<bool, !std::is_convertible<TType,int>{}
                                                            && std::is_enum<TType>{}>;

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor
    /// @tparam TKeyType of this dictionary level.
    /// @param rKey Key of this dictionary level.
    /// @param rKeyName Name of the Key. Used to print the Key name.
    /// @param IsList Optional parameter to declare this dictionary a list. Default: false.
    template<typename TKeyType>
    Dictionary(TKeyType rKey, std::string rKeyName, bool IsList=false) :
        mKey(rKey), mKeyName(rKeyName), mKeyTypeName(GetTypeName<TKeyType>()),
        mChildrenDummyKey(std::monostate{}), mChildrenKeyTypeName(""), mIsList(IsList)
    {
    }

    /// Destructor
    virtual ~Dictionary() = default;
    /// Copy Constructor
    Dictionary(const Dictionary& rDict) = default;
    /// Assignement operator
    Dictionary& operator=(const Dictionary& rDict) = default;
    /// Move constructor
    Dictionary(Dictionary&& rDict) = default;
    /// Move assignement operator
    Dictionary& operator=(Dictionary&& rDict) = default;

    ///@}
    ///@name Operations
    ///@{

    ////////////////
    /// AddItems ///
    ////////////////

    /// @brief Adds an empty subdictionary.
    /// @tparam TKeyType
    /// @param rKey
    /// @param rKeyName
    /// @return Dictionary&
    template<typename TKeyType>
    Dictionary& AddEmptySubDictionary(TKeyType rKey, std::string rKeyName) {
        static_assert( is_scoped_int_enum<TKeyType>::value );
        // Check if key index fits the current dictionary size.
        const int index = static_cast<int>(rKey) - DictStarts::start_subdicts;
        QuESo_ERROR_IF( static_cast<int>(mSubDictionaries.size()) != index)
            << "Dictionaries are not added in correct order.\n";

        if( mChildrenDummyKey.index() == 0 ) { // Dummy key not set -> std::(monostate: index()==0).
            mChildrenDummyKey = rKey; // Store key to be able to test the type later.
            mChildrenKeyTypeName = GetTypeName<TKeyType>();
        }
        else {
            if( !std::get_if<TKeyType>(&mChildrenDummyKey) ) {
                QuESo_ERROR << "The keys on the same level must be of same type. Given Key type: '" << GetTypeName<TKeyType>()
                    << ". Stored Key type on this level: '" << mChildrenKeyTypeName << "'.\n";
            }
        }

        // Add new dictionary
        mSubDictionaries.push_back( Dictionary(rKey, rKeyName) );
        return mSubDictionaries.back();
    }

    ///@brief Adds an item to the list.
    ///@tparam TKeyType
    ///@tparam TValueTypes
    ///@param NewEntries
    ///@return Dictionary&
    template<typename TKeyType, typename... TValueTypes>
    Dictionary& AddListItem(std::tuple<std::tuple<TKeyType, std::string, TValueTypes, bool>...> NewEntries) {
        static_assert( is_scoped_int_enum<TKeyType>::value );
        QuESo_ERROR_IF(!mIsList) << "Trying to add a ListItem to an object, which is not a list.\n";
        // A list item is again a list.
        mListsOfDicts.push_back( Dictionary(std::monostate{}, "") );
        auto& r_new_list = mListsOfDicts.back();
        r_new_list.AddValues(NewEntries);
        return r_new_list;
    }

    /// @brief Adds an empty subdictionary.
    /// @tparam TKeyType
    /// @param rKey
    /// @param rKeyName
    /// @return Dictionary&
    template<typename TKeyType>
    Dictionary& AddEmptyList(TKeyType rKey, std::string rKeyName) {
        static_assert( is_scoped_int_enum<TKeyType>::value );
        // Check if key index fits the current dictionary size.
        const int index = static_cast<int>(rKey) - DictStarts::start_lists;
        QuESo_ERROR_IF(static_cast<int>(mListsOfDicts.size()) != index)
            << "Lists are not added in correct order.\n";

        if( mChildrenDummyKey.index() == 0 ) { // Dummy key not set -> std::(monostate: index()==0).
            mChildrenDummyKey = rKey; // Store key to be able to test the type later.
            mChildrenKeyTypeName = GetTypeName<TKeyType>();
        }
        else {
            if( !std::get_if<TKeyType>(&mChildrenDummyKey) ) {
                QuESo_ERROR << "The keys on the same level must be of same type. Given Key type: '" << GetTypeName<TKeyType>()
                    << ". Stored Key type on this level: '" << mChildrenKeyTypeName << "'.\n";
            }
        }

        // Add new list
        constexpr bool is_list = true;
        mListsOfDicts.push_back( Dictionary(rKey, rKeyName, is_list) );
        return mListsOfDicts.back();
    }

    /// @brief Adds values to the current level of the dictionary. Should only be called in derived class to define the default dictionary.
    ///        Layout: std::tuple( std::tuple(Enum Key, string KeyName, ValueType Value, bool IsSet ) )
    ///        If IsSet=false, 'Value' is only used to deduce the associated type.
    /// @tparam TKeyType of keys. All keys within one level must have the same type (same enum).
    /// @tparam ...TValueTypes Parameter pack for types.
    /// @param NewEntries
    template<typename TKeyType, typename... TValueTypes>
    void AddValues(std::tuple<std::tuple<TKeyType, std::string, TValueTypes, bool>...> NewEntries) {
        static_assert( is_scoped_int_enum<TKeyType>::value );
        if( mChildrenDummyKey.index() == 0 ) { // Dummy key not set -> std::(monostate: index()==0).
            // We store one key, to be able to test for the Key Type.
            mChildrenDummyKey = std::get<0>(std::get<0>(NewEntries));
            mChildrenKeyTypeName = GetTypeName<TKeyType>();
        }
        UnpackTuple(NewEntries);
    }

    ////////////////
    /// GetItems ///
    ////////////////

    /// @brief Returns ref to subdictionary corresponding to given key provided as Enum (const version).
    /// @tparam TKeyType.
    /// @param QueryKey
    /// @return const Dictionary&
    /// @note only throws in DEBUG mode.
    template<typename TKeyType,
             typename = std::enable_if_t<is_scoped_int_enum<TKeyType>::value>>
    const Dictionary& operator [](TKeyType QueryKey) const noexcept(NOTDEBUG) {
        QuESo_ASSERT( std::get_if<TKeyType>(&mChildrenDummyKey) != 0, "Given Key type does not match stored Key type.\n");
        const int index = static_cast<int>(QueryKey) - DictStarts::start_subdicts;
        return mSubDictionaries[index];
    }

    /// @brief Returns ref to subdictionary corresponding to given key provided as Enum (non const version).
    /// @tparam TKeyType.
    /// @param QueryKey
    /// @return Dictionary&
    /// @note only throws in DEBUG mode.
    template<typename TKeyType,
             typename = std::enable_if_t<is_scoped_int_enum<TKeyType>::value>>
    Dictionary& operator [](TKeyType QueryKey) noexcept(NOTDEBUG) {
        return const_cast<Dictionary&>(static_cast<const Dictionary&>(*this)[QueryKey]);
    }

    /// @brief Returns ref to subdictionary corresponding to given key provided as std::string.
    ///        This function should only be used in Python. The Version: [Enum Key] is much faster.
    /// @param rQueryKeyName
    /// @return Dictionary&
    Dictionary& operator[] (const std::string& rQueryKeyName) {
        // Lets search for the key name.
        const auto it = std::find_if(mSubDictionaries.begin(), mSubDictionaries.end(),
            [&rQueryKeyName](const Dictionary& rDict)
                { return rQueryKeyName == rDict.GetKeyName();});
        if( it == mSubDictionaries.end() ){
            QuESo_ERROR << "Given Key (" << rQueryKeyName << ") was not found.\n";
        }
        IndexType index = it - mSubDictionaries.begin();
        return mSubDictionaries[index];
    }

    /// @brief Get the List to the given Key provided as Enum (const version).
    /// @tparam TKeyType
    /// @param QueryKey
    /// @return const std::vector<Dictionary>&
    /// @note only throws in DEBUG mode.
    template<typename TKeyType>
    const std::vector<Dictionary>& GetList(TKeyType QueryKey,
                                           std::enable_if_t<is_scoped_int_enum<TKeyType>::value>* = nullptr ) const noexcept(NOTDEBUG) {
        QuESo_ASSERT( std::get_if<TKeyType>(&mChildrenDummyKey) != 0, "Given Key Type does not match stored Key type.\n" );
        const int index = static_cast<int>(QueryKey) - DictStarts::start_lists;
        return mListsOfDicts[index].GetList();
    }

    /// @brief Get the List to the given Key provided as Enum (non const version).
    /// @tparam TKeyType
    /// @param QueryKey
    /// @return std::vector<Dictionary>&
    /// @note only throws in DEBUG mode.
    template<typename TKeyType>
    std::vector<Dictionary>& GetList(TKeyType QueryKey,
                                     std::enable_if_t<is_scoped_int_enum<TKeyType>::value>* = nullptr ) noexcept(NOTDEBUG) {
        return const_cast<std::vector<Dictionary>&>(static_cast<const Dictionary&>(*this).GetList(QueryKey));
    }

    /// @brief Get the List to the given Key provided as std::string.
    ///        This function should only be used in Python. The Version: [Enum Key] is much faster.
    /// @param rQueryKeyName
    /// @return std::vector<Dictionary>&
    std::vector<Dictionary>& GetList(const std::string& rQueryKeyName) {
        // Lets search for the key name.
        const auto it = std::find_if(mListsOfDicts.begin(), mListsOfDicts.end(),
            [&rQueryKeyName](const Dictionary& rDict)
                { return rQueryKeyName == rDict.GetKeyName();});
        if( it == mListsOfDicts.end() ){
            QuESo_ERROR << "Given Key (" << rQueryKeyName << ") was not found.\n";
        }
        IndexType index = it - mListsOfDicts.begin();
        return mListsOfDicts[index].GetList();
    }

    ///@brief Get the List Object to a given key. In contrast to GetList(), this functions returns the entire Dictionary,
    ///       which holds the key and the underlying list.
    ///@tparam TKeyType
    ///@param QueryKey
    ///@return Dictionary&
    ///@see GetList()
    template<typename TKeyType>
    Dictionary& GetListObject(TKeyType QueryKey) {
        static_assert( is_scoped_int_enum<TKeyType>::value );
        if ( std::get_if<TKeyType>(&mChildrenDummyKey) ){
            const IndexType index = static_cast<IndexType>(QueryKey) - DictStarts::start_lists;
            return mListsOfDicts[index];
        }
        QuESo_ERROR << "Given Key type (" << GetTypeName<TKeyType>() << ") does not match stored Key type (" << mChildrenKeyTypeName << ").\n";
    }

    /// @brief Returns the Value to the given Key provided as Enum.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param QueryKey (Enum)
    /// @return const TValueType&
    /// @note Only throws in DEBUG mode.
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValue(TKeyType QueryKey,
                               std::enable_if_t<is_scoped_int_enum<TKeyType>::value>* = nullptr ) const noexcept(NOTDEBUG){
        QuESo_ASSERT( std::get_if<TKeyType>(&mChildrenDummyKey) != 0, "Given Key Type does not match stored Key type.\n" );
        const int index = static_cast<int>(QueryKey) - DictStarts::start_values;
        return mData[index].template GetValueNoCheck<TValueType>();
    }

    /// @brief Returns the Value to the given Key provided as std::string.
    ///        This function should only be used in Python. The Version: [Enum Key] is much faster.
    /// @tparam TValueType
    /// @param QueryKey (std::string)
    /// @return const TValueType&
    template<typename TValueType>
    const TValueType& GetValue(const std::string& rQueryKeyName) const {
        IndexType index = FindIndexOfKeyValuePairByString(rQueryKeyName);
        return mData[index].template GetValue<TValueType>();
    }

    ////////////////
    /// SetItems ///
    ////////////////

    /// @brief Sets Value to given Key provided as Enum.
    /// @tparam TKeyType
    /// @tparam TValueType
    /// @param QueryKey (Enum)
    /// @param NewValue
    /// @see Only throws in DEBUG mode.
    template<typename TKeyType, typename TValueType>
    void SetValue(TKeyType QueryKey, TValueType NewValue,
                  std::enable_if_t<is_scoped_int_enum<TKeyType>::value>* = nullptr ) noexcept(NOTDEBUG) {
        QuESo_ASSERT( std::get_if<TKeyType>(&mChildrenDummyKey) != 0, "Given Key Type does not match stored Key type.\n" );
        const IndexType index = static_cast<IndexType>(QueryKey) - DictStarts::start_values;
        if constexpr(std::is_same_v<TValueType, int>
                  || std::is_same_v<TValueType, unsigned int>
                  || std::is_same_v<TValueType, unsigned long>) { // Accept different integer types.
            if constexpr (std::is_same_v<TValueType, int>) {
                QuESo_ASSERT(NewValue >= 0, "Value must be non-negative.\n");
            }
            mData[index].SetValueNoCheck(static_cast<IndexType>(NewValue));
        } else {
            mData[index].SetValueNoCheck(NewValue);
        }
    }

    /// @brief Sets Value to given Key (given as string).
    /// @tparam TValueType
    /// @param rQueryKeyName
    /// @param NewValue
    /// @see SetValueWithAmbiguousType() <- allows to cast ambiguous types, e.g., casts 0 -> 0.0, if possible.
    /// @note Should only be used in Python. There is also a version that takes an actual KeyType instead of a string.
    template<typename TValueType>
    void SetValue(const std::string& rQueryKeyName, TValueType NewValue) {
        IndexType index = FindIndexOfKeyValuePairByString(rQueryKeyName);
        if constexpr(std::is_same_v<TValueType, int>
                || std::is_same_v<TValueType, unsigned int>
                || std::is_same_v<TValueType, unsigned long>) { // Accept different integer types.
            if constexpr (std::is_same_v<TValueType, int>) {
                QuESo_ERROR_IF(NewValue < 0) << "Value must be non-negative.\n";
            }
            mData[index].SetValue(static_cast<IndexType>(NewValue));
        } else {
            mData[index].SetValue(NewValue);
        }
    }

    /// @brief Sets Value to given Key (given as string). Casts 0 -> 0.0 and [0, 0, 0] -> [0.0, 0.0, 0.0] if possible.
    /// @tparam TValueType
    /// @param rQueryKeyName
    /// @param NewValue
    /// @param SetValue() <- Version without type casting.
    /// @note Should only be used in Python. There is also a version that takes an actual KeyType instead of a string.
    template<typename TValueType>
    void SetValueWithAmbiguousType(const std::string& rQueryKeyName, TValueType NewValue) {
        IndexType index = FindIndexOfKeyValuePairByString(rQueryKeyName);
        if constexpr(std::is_same_v<TValueType, int>
                || std::is_same_v<TValueType, unsigned int>
                || std::is_same_v<TValueType, unsigned long>) { // Accept different integer types.
            if constexpr (std::is_same_v<TValueType, int>) {
                QuESo_ERROR_IF(NewValue < 0) << "Value must be non-negative.\n";
            }
            mData[index].SetValueWithAmbiguousType(static_cast<IndexType>(NewValue));
        } else {
            mData[index].SetValueWithAmbiguousType(NewValue);
        }
    }

    ///////////////////////////
    /// AdditionalFunctions ///
    ///////////////////////////

    /// @brief Returns true if Value to given Key is set.
    /// @tparam TKeyType.
    /// @param QueryKey
    /// @return bool
    /// @note Only throws in DEBUG mode.
    template<typename TKeyType>
    bool IsSet(TKeyType QueryKey,
               std::enable_if_t<is_scoped_int_enum<TKeyType>::value>* = nullptr) const noexcept(NOTDEBUG) {
               static_assert( is_scoped_int_enum<TKeyType>::value );
        QuESo_ASSERT( std::get_if<TKeyType>(&mChildrenDummyKey) != 0, "Given Key Type does not match stored Key type.\n" );
        const int index = static_cast<int>(QueryKey) - DictStarts::start_values;
        return mData[index].IsSet();
    }

    /// @brief Returns true if Value to given Key is set.
    /// @param rQueryKeyName
    /// @return bool
    bool IsSet(const std::string& rQueryKeyName) const {
        IndexType index = FindIndexOfKeyValuePairByString(rQueryKeyName);
        return mData[index].IsSet();
    }

    /// @brief Throws error if any stored value on this level is not set.
    void CheckIfValuesAreSet() const {
        for( const auto& r_key_val_pair : mData  ){
            QuESo_ERROR_IF( !r_key_val_pair.IsSet() ) << "Variable: " << r_key_val_pair.GetKeyName() << " is not set.\n";
        }
    }

    /// @brief Prints this dictionary in JSON format.
    /// @param rOStream
    void PrintInfo(std::ostream& rOStream) const {
        std::string indent = "";
        PrintInfo(rOStream, indent);
    }

private:

    ///@}
    ///@name Private operations.
    ///@{

    /// @brief Returns Key of current dictionary level
    /// @return  const TVariantKey&
    const TVariantKey& GetKey() const noexcept {
        return mKey;
    }

    /// @brief Returns Key of current dictionary level
    /// @return  const TVariantKey&
    const std::string& GetKeyName() const noexcept {
        return mKeyName;
    }

    /// @brief Returns name of Key type.
    /// @return const std::string&
    const std::string& GetKeyTypeName() const noexcept {
        return mKeyTypeName;
    }

    ///@brief Get the List object (const version).
    ///@return const std::vector<Dictionary>&
    const std::vector<Dictionary>& GetList() const noexcept {
        return mListsOfDicts;
    }

    ///@brief Get the List object (non-const version).
    ///@return const std::vector<Dictionary>&
    std::vector<Dictionary>& GetList() noexcept {
        return mListsOfDicts;
    }

    /// Returns type name of TType as char*
    template<typename TType>
    const char* GetTypeName() const {
        // Names are not really pretty, but abi::__cxa__demangle only works for gcc.
        return typeid(TType).name();
    }

    ///@brief Returns index of KeyValuePair to given Key provided as std::string.
    ///@param rQueryKeyName Key.
    ///@return IndexType
    IndexType FindIndexOfKeyValuePairByString(const std::string& rQueryKeyName) const {
        // Lets search for the key name.
        const auto it = std::find_if(mData.begin(), mData.end(),
            [&rQueryKeyName](const KeyValuePairType& rValue)
                { return rQueryKeyName == rValue.GetKeyName();});
        if( it != mData.end() ){
            IndexType index = it - mData.begin();
            return index;
        }
        QuESo_ERROR << "Given Key (" << rQueryKeyName << ") was not found.\n";
    }

    /// @brief Prints this dictionary in JSON format.
    /// @param rOStream
    /// @param rIndent This variable is used for recursive function calls.
    /// @param IsLastItem This variable is used for recursive function calls.
    void PrintInfo(std::ostream& rOStream, std::string& rIndent, bool IsLastItem=false) const {
        /// Add proper header (with old indent).
        if( rIndent == "" ){
            rOStream << "{\n"; // Start root
        } else if (mKey.index() == 0) { // Key is in std::monostate
            rOStream << rIndent << "{\n"; // Start list item
        } else if (mIsList) {
            rOStream << rIndent << '\"' << mKeyName << "\" : [\n"; // Start List
        } else {
            rOStream << rIndent << '\"' << mKeyName << "\" : {\n"; // Start dictionary
        }

        std::string new_indent = rIndent + "\t";
        /// Write KeyValuePairs with new indent.
        for( IndexType i = 0; i < mData.size(); ++i){
            const auto& r_key_value_pair = mData[i];
            rOStream << new_indent << "\"" << r_key_value_pair.GetKeyName() << "\" : ";
            r_key_value_pair.PrintValue(rOStream);
            if( i == mData.size() - 1 && mSubDictionaries.size() == 0 && mListsOfDicts.size() == 0){
                rOStream << '\n';
            } else {
                rOStream << ",\n";
            }
        }

        /// Write mSubDictionaries with new indent.
        for( IndexType i = 0; i < mSubDictionaries.size(); ++i){
            if( i == mSubDictionaries.size() - 1 && mListsOfDicts.size() == 0) {
                mSubDictionaries[i].PrintInfo(rOStream, new_indent, true);
            } else {
                mSubDictionaries[i].PrintInfo(rOStream, new_indent, false);
            }
        }

        /// Write mListsOfDicts with new indent.
        for( IndexType i = 0; i < mListsOfDicts.size(); ++i){
            if( i == mListsOfDicts.size() - 1) {
                mListsOfDicts[i].PrintInfo(rOStream, new_indent, true);
            } else {
                mListsOfDicts[i].PrintInfo(rOStream, new_indent, false);
            }
        }

        /// Close with old indent.
        if( rIndent == "" || IsLastItem ){ // End list or dict item.
            rOStream << ((mIsList) ? (rIndent + "]") : (rIndent + "}")) << "\n"; // Needs to be encapsulated!
        } else {
            rOStream << ((mIsList) ? (rIndent + "],") : (rIndent + "},")) << "\n";
        }
    }

    ///@}
    ///@name Private struct and function definitions to extract tuples.
    ///@{

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
        if(  (static_cast<int>(std::get<0>(r_tuple)) - DictStarts::start_values) != I ){
            QuESo_ERROR << "Tuples are not given in correct order. \n";
        }

        mData.push_back(KeyValuePairType(TVariantKey(std::get<0>(r_tuple)), std::get<1>(r_tuple), std::get<2>(r_tuple), std::get<3>(r_tuple)));
        UnpackTuple<I+1>(Tuple); // Recursive call.
    }

    ///@}
    ///@name Private member variables.
    ///@{

    /// Key related members. (Key of this dictionary)
    TVariantKey mKey;
    std::string mKeyName;
    std::string mKeyTypeName;

    /// Members related to the key related to all possible children.
    // Note that mData, mSubDictionaries, mListsOfDicts use the same key types.
    TVariantKey mChildrenDummyKey; // This dummy key allows to test if any given Key type matches the stored Key types.
    std::string mChildrenKeyTypeName;

    /// Vector of KeyValuePairs.
    std::vector<KeyValuePairType> mData;

    /// Vector of Subdictionaries.
    std::vector<Dictionary<TEnumKeys...>> mSubDictionaries;

    /// Members related to Lists.
    bool mIsList;
    std::vector<Dictionary<TEnumKeys...>> mListsOfDicts;

    ///@}
};

/// Output stream functions
template<typename... TEnumKeys>
inline std::ostream& operator<< (std::ostream& rOStream, const Dictionary<TEnumKeys ...>& rThis) {
    rThis.PrintInfo(rOStream);
    return rOStream;
}
///@} // End QuESo Classes
} // End queso namespace.

#endif // End DICTIONARY_INCLUDE_HPP