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
#include <cassert>
#include <vector>
#include <algorithm>
#include <type_traits>

/// Project includes
#include "queso/containers/data_set.hpp"
#include "queso/includes/define.hpp"

namespace queso {

/// Forward declaration
template<typename DictionaryType>
class DictionaryStringAccess;

///@name QuESo Classes
///@{


/// @class  Dictionary.
/// @author Manuel Messmer
/// @tparam TKeySetValuesTypeTag: Defines the possible types that can be stored in Dictionary, e.g., double, IndexType, etc.
/// @brief Dictionary for fast data access. 'TKeySetValuesTypeTag' must be registered in 'register_keys.hpp'.
///        The DataSet allows to access data quickly by ConstexprKeys -> 'register_keys.hpp', also see 'keys.hpp'.
/// @details Dictionary can contain DataSets (KeyValuePairs), SubDictionaries, and Lists. For each 'section', a unique KeySet must
///        be registered in 'register_keys.hpp'. However, the keys for all three groups are made available within one namespace.
///        An example KeySet with ConstexprKeys looks like this:
///        - KeySetName::value_one -> Used in, e.g., GetValue().
///        - KeySetName::value_two -> Used in, e.g., GetValue().
///        - KeySetName::sub_dict_one -> Used in, e.g., operator[] (return SubDictionary).
///        - KeySetName::sub_dict_two -> Used in, e.g., operator[] (return SubDictionary).
///        - KeySetName::list_one -> Used in, e.g., GetList().
///        - KeySetName::list_two -> Used in, e.g., GetList().
///        Each key contains a TypeTag, e.g., ValueTypeTag (template param), SubDictionaryTypeTag, ListTypeTag.
///        With these type tags, we can make sure that each function, e.g., GetValue() is actually called with a correct key:
///
///        Besides accessing the data via ConstexprKeys, the DataSet also provides methods to access the data via std::strings.
///        Given the KeyName a DynamicKey is deduced from KeySetInfo (which is passed to the constructor). The respective KeySetInfo
///        is also registered in 'register_keys.hpp'. The main idea is to access the data via ConstexprKeys in C++.
///        The std::string access is meant to be used in Python only. Therefore, the respective functions are private and made
///        available via DictionaryStringAccess. This is just a safety guard that these functions are not accidentally used in C++.
/// @see 'register_keys.hpp' and 'keys.hpp'.
template<typename TKeySetValuesTypeTag>
class Dictionary
{
public:
    ///@name Type definitions
    ///@{

    using KeySetValuesTypeTag = TKeySetValuesTypeTag;
    using DataSetType = DataSet<KeySetValuesTypeTag>;
    using ListType = std::vector<Unique<Dictionary>>;

    /// Helper to check if TType is a key.
    template<typename TType>
    static constexpr bool is_key_v = KeySetValuesTypeTag::template is_key_v<TType>;

    /// Type trait helpers.
    template<typename TType>
    static constexpr bool has_value_keys_v = !std::is_same_v<typename TType::ValueKeySetInfoType, std::false_type>;

    template<typename TType>
    static constexpr bool has_subdict_keys_v = !std::is_same_v<typename TType::SubDictKeySetInfoType, std::false_type>;

    template<typename TType>
    static constexpr bool has_list_keys_v = !std::is_same_v<typename TType::ListKeySetInfoType, std::false_type>;

    /// Helper type tag to use as constructor template.
    template<typename TValueKeySetInfoType,
             typename TSubDictKeySetInfoType,
             typename TListKeySetInfoType>
    struct KeySetInfosTypeTag {
        using ValueKeySetInfoType = TValueKeySetInfoType;
        using SubDictKeySetInfoType = TSubDictKeySetInfoType;
        using ListKeySetInfoType = TListKeySetInfoType;
    };

    /// Helper to create an empty KeySet.
    using EmptyKeySetType = std::false_type;

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor. Dictionary takes one, two, or three KeySetInfoType's. These can be passed using
    ///        KeySetInfosTypeTag. Each KeySetInfoType contains the keys to access Values, SubDictionaries and Lists.
    ///        Each group can also be empty (use EmptyKeySetType for the respective KeySetInfoType).
    ///        Example use:
    ///            using DictionaryType = Dictionary<key::ValuesTypeTag>;
    ///            DictionaryType test_dict( DictionaryType::KeySetInfosTypeTag<
    ///                key::detail::KeySetNameValuesTypeTagKeySetInfo,
    ///                key::detail::KeySetNameSubDictTypeTagKeySetInfo,
    ///                DictionaryType::EmptyKeySetType>{} ); // <- List section would be empty.
    /// @tparam TKeySetInfosTypeTag
    /// @param SetInfoTypeTag
    template<typename TKeySetInfosTypeTag>
    Dictionary( TKeySetInfosTypeTag ) {
        using Tag = TKeySetInfosTypeTag;

        if constexpr (has_value_keys_v<Tag>) { // Initialize DataSet.
            using KetSetInfoType = typename Tag::ValueKeySetInfoType;
            using DataSetTypeTag = typename DataSetType::template KeySetInfoTypeTag<KetSetInfoType>;
            mpDataSet = MakeUnique<DataSetType>( DataSetTypeTag{} );
        }
        if constexpr ( has_subdict_keys_v<Tag> ) { // Initialize SubDictionaries.
            using KetSetInfoType = typename Tag::SubDictKeySetInfoType;
            mpSubDictKeySetInfo = MakeUnique<KetSetInfoType>();
            mSubDictionaries.resize(mpSubDictKeySetInfo->GetNumberOfKeys());
        }
        if constexpr ( has_list_keys_v<Tag> ) { // Initialize Lists.
            using KetSetInfoType = typename Tag::ListKeySetInfoType;
            mpListKeySetInfo = MakeUnique<KetSetInfoType>();
            mListsOfDicts.resize(mpListKeySetInfo->GetNumberOfKeys());
        }
    }

    /// Destructor
    ~Dictionary() = default;
    /// Copy Constructor
    Dictionary(const Dictionary& rDict) = delete;
    /// Assignement operator
    Dictionary& operator=(const Dictionary& rDict) = delete;
    /// Move constructor
    Dictionary(Dictionary&& rDict) = default;
    /// Move assignement operator
    Dictionary& operator=(Dictionary&& rDict) = default;

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
    void SetValue(const TKeyType& rQueryKey, const TValueType& rNewValue) noexcept(NOTDEBUG) {
        QuESo_ASSERT( mpDataSet != nullptr, "This dictionary has an empty data set.\n");
        mpDataSet->SetValue(rQueryKey, rNewValue);
    }

    /// @brief Returns the value associated with the given Key.
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return const TValueType&
    /// @note In release mode it throws if the given key is not set. In debug mode it addionally throws if TKeyType is of wrong type.
    /// @see GetValueFast(): Never throws in release.
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValue(const TKeyType& rQueryKey) const {
        QuESo_ASSERT( mpDataSet != nullptr, "This dictionary has an empty data set.\n");
        return mpDataSet->template GetValue<TValueType>(rQueryKey);
    }

    /// @brief Returns the value associated with the given Key (fast version, does not throw).
    /// @tparam TValueType
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return const TValueType&
    /// @note Only throws in debug mode. In realese it does not check if the value is set and if TKeyType is of correct type.
    /// @see GetValue(): Throws if value is not set.
    template<typename TValueType, typename TKeyType>
    const TValueType& GetValueFast(const TKeyType& rQueryKey) const noexcept(NOTDEBUG) {
        QuESo_ASSERT( mpDataSet != nullptr, "This dictionary has an empty data set.\n");
        return mpDataSet->template GetValueFast<TValueType>(rQueryKey);
    }

    /// @brief Returns true if the value associated with the given Key is set.
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return bool
    /// @note Only throws in debug mode.
    template<typename TKeyType>
    bool IsSet(const TKeyType& rQueryKey) const noexcept(NOTDEBUG) {
        QuESo_ASSERT( mpDataSet != nullptr, "This dictionary has an empty data set.\n");
        return mpDataSet->IsSet(rQueryKey);
    }

    /// @brief Sets a SubDictionary associated with the given key. Dictionary takes unique ownership of the SubDictionary.
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @param rNewDictionary
    /// @note Only throws in debug mode.
    template<typename TKeyType>
    void SetSubDictionary(const TKeyType& rQueryKey, Unique<Dictionary>&& rNewDictionary,
                           std::enable_if_t<is_key_v<TKeyType>>* = nullptr  ) noexcept(NOTDEBUG) {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeySetToWhat, key::SubDictTypeTag>::value );
        static_assert( std::is_same<typename TKeyType::KeyToWhat, key::SubDictTypeTag>::value );

        QuESo_ASSERT( mpSubDictKeySetInfo != nullptr, "This dictionary does not contain any sub-dictionaries.\n");
        QuESo_ASSERT( mpSubDictKeySetInfo->IsSameKeySet(rQueryKey.KeySetInfoTypeIndex()),
            "Given Key: '" + std::string(rQueryKey.Name()) + "' is of wrong type.\n" );

        mSubDictionaries[rQueryKey.Index()] = std::move(rNewDictionary);
    }

    /// @brief Returns SubDictionary associated with the given Key (const version).
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return const Dictionary&
    /// @note Only throws in debug mode.
    template<typename TKeyType,
             typename = std::enable_if_t<is_key_v<TKeyType>>>
    const Dictionary& operator[](const TKeyType& rQueryKey) const noexcept(NOTDEBUG) {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeySetToWhat, key::SubDictTypeTag>::value );
        static_assert( std::is_same<typename TKeyType::KeyToWhat, key::SubDictTypeTag>::value );

        QuESo_ASSERT( mpSubDictKeySetInfo != nullptr, "This dictionary does not contain any sub-dictionaries.\n");
        QuESo_ASSERT( mpSubDictKeySetInfo->IsSameKeySet(rQueryKey.KeySetInfoTypeIndex()),
            "Given Key: '" + std::string(rQueryKey.Name()) + "' is of wrong type.\n");

        QuESo_ASSERT( rQueryKey.Index() < mSubDictionaries.size(), "It seems like, we forgot to resize 'mSubDictionaries'.\n");
        QuESo_ASSERT( mSubDictionaries[rQueryKey.Index()] != nullptr,
            "Dictionary associated with Key: " + std::string( rQueryKey.Name() ) + " is not set.\n");

        return *(mSubDictionaries[rQueryKey.Index()]);
    }

    /// @brief Returns SubDictionary associated with the given Key (non-const version).
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return Dictionary&
    /// @note Only throws in debug mode.
    template<typename TKeyType,
             typename = std::enable_if_t<is_key_v<TKeyType>>>
    Dictionary& operator[](const TKeyType& rQueryKey) noexcept(NOTDEBUG) {
        return const_cast<Dictionary&>( static_cast<const Dictionary&>(*this)[rQueryKey] );
    }

    /// @brief Returns List associated with the given Key (const version).
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return const ListType&
    /// @note Only throws in debug mode.
    template<typename TKeyType>
    const ListType& GetList(const TKeyType& rQueryKey,
                            std::enable_if_t<is_key_v<TKeyType>>* = nullptr  ) const noexcept(NOTDEBUG) {

        static_assert( std::is_same<typename TKeyType::KeySetInfoType::KeySetToWhat, key::ListTypeTag>::value );
        static_assert( std::is_same<typename TKeyType::KeyToWhat, key::ListTypeTag>::value );

        QuESo_ASSERT( mpListKeySetInfo != nullptr, "This dictionary does not contain any lists.\n");

        QuESo_ASSERT( mpListKeySetInfo->IsSameKeySet(rQueryKey.KeySetInfoTypeIndex()),
            "Given Key: '" + std::string(rQueryKey.Name()) + "' is of wrong type.\n" );

        QuESo_ASSERT( rQueryKey.Index() < mListsOfDicts.size(), "It seems like, we forgot to resize 'mListsOfDicts'.\n");

        return mListsOfDicts[rQueryKey.Index()];
    }

    /// @brief Returns List associated with the given Key (non-const version).
    /// @tparam TKeyType
    /// @param rQueryKey
    /// @return const ListType&
    /// @note Only throws in debug mode.
    template<typename TKeyType>
    ListType& GetList(const TKeyType& rQueryKey,
                      std::enable_if_t<is_key_v<TKeyType>>* = nullptr  ) noexcept(NOTDEBUG) {
        return const_cast<ListType&>(static_cast<const Dictionary&>(*this).GetList(rQueryKey));
    }

private:
    ///@}
    ///@name Type definitions
    ///@{

    using DataStringAccess = DataSetStringAccess<DataSetType>;

    ///@}
    ///@name Private operations
    ///@{

    /// Friend to access private member: PrintInfo.
    template<typename Tag>
    friend std::ostream& operator<< (std::ostream& rOStream, const Dictionary<Tag>& rDictionary);

    /// @brief Prints this dictionary (recursively) in JSON format.
    /// @param rOStream
    /// @param rIndent This variable is used for recursive function calls.
    void PrintInfo(std::ostream& rOStream, const std::string& rIndent) const {
        std::string new_indent = rIndent + "    ";

        const IndexType n_subdicts = mSubDictionaries.size();
        const IndexType n_lists = mListsOfDicts.size();
        const IndexType n_items = n_subdicts + n_lists;

        /// Write DataSet to rOstream.
        if( mpDataSet ) {
            const auto& r_data_key_set_info = mpDataSet->GetKeySetInfo();
            const auto& r_map = r_data_key_set_info.GetStringKeyMap();
            auto it = r_map.begin();
            while (it != r_map.end()) {
                const auto& [name, p_key] = *it;
                rOStream << new_indent << "\"" << name << "\" : ";
                mpDataSet->PrintValue(p_key, rOStream);
                rOStream << ((++it == r_map.end() && n_items == 0) ? "\n" : ",\n");
            }
        }

        /// Write SubDictionaties (recursively) to rOStream.
        if( mpSubDictKeySetInfo ) {
            const auto& r_map = mpSubDictKeySetInfo->GetStringKeyMap();
            auto it = r_map.begin();
            while (it != r_map.end()) {
                const auto& [r_name, p_key] = *it;
                const auto& p_sub_dict = mSubDictionaries[p_key->Index()];
                rOStream << new_indent << '\"' << r_name << "\" : {\n"; // Start dictionary.
                if( p_sub_dict ) {
                    p_sub_dict->PrintInfo(rOStream, new_indent);
                }
                if( ++it == r_map.end() && n_lists == 0 ) { // End dictionary.
                   rOStream << new_indent << "}\n"; // Is last item.
                } else {
                    rOStream << new_indent << "},\n"; // Is not last item.
                }
            }
        }

        /// Write Lists (recursively) to rOStream.
        if( mpListKeySetInfo ) {
            const auto& r_map = mpListKeySetInfo->GetStringKeyMap();
            auto it = r_map.begin();
            while (it != r_map.end()) {
                const auto& [r_name, p_key] = *it;
                const auto& r_list = mListsOfDicts[p_key->Index()];
                rOStream << new_indent << '\"' << r_name << "\" : [\n"; // Start List.
                std::string inner_list_indent = new_indent + "    ";
                auto it_inner = r_list.begin();
                while ( it_inner != r_list.end() ){
                    rOStream << inner_list_indent << "{\n"; // Start list item.
                    (*it_inner)->PrintInfo(rOStream, inner_list_indent);
                    if( ++it_inner ==  r_list.end() ){ // End list item.
                        rOStream << inner_list_indent << "}\n"; // Is last item.
                    } else {
                        rOStream << inner_list_indent << "},\n"; // Is not last item.
                    }
                }
                if( ++it == r_map.end() ) { // End list.
                    rOStream << (new_indent + "]\n"); // Is last item.
                } else {
                    rOStream << (new_indent + "],\n"); // Is not last item.
                }
            }
        }
    }

    /// The following operations allow to access the data via KeyNames (std::string).
    /// These functions should only be used in Python. Therefore, they are made private here
    /// and can only be used via the friend class DictionaryStringAccess.
    friend class DictionaryStringAccess<Dictionary>;

    /// @brief Sets a SubDictionary associated with the given KeyName. Dictionary takes unique ownership of the SubDictionary.
    /// @param rQueryKeyName (std::string).
    /// @param rNewDictionary
    void SetSubDictionary(const std::string& rQueryKeyName, Unique<Dictionary>&& rNewDictionary) {

        QuESo_ERROR_IF( !mpSubDictKeySetInfo )  << "This dictionary does not contain any sub-dictionaries.\n";

        const auto p_key = mpSubDictKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName
                               << "'. Possible names are: " << mpSubDictKeySetInfo->GetAllKeyNames() << "\n.";

        QuESo_ASSERT( p_key->Index() < mSubDictionaries.size(), "It seems like, we forgot to resize 'mSubDictionaries'.\n");

        mSubDictionaries[p_key->Index()] = std::move(rNewDictionary);
    }

    /// @brief Returns SubDictionary associated with the given KeyName.
    /// @param rQueryKeyName (std::string).
    /// @return Dictionary&
    Dictionary& operator[] (const std::string& rQueryKeyName) {

        QuESo_ERROR_IF( !mpSubDictKeySetInfo )  << "This dictionary does not contain any sub-dictionaries.\n";

        const auto p_key = mpSubDictKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName
                               << "'. Possible names are: " << mpSubDictKeySetInfo->GetAllKeyNames() << "\n.";

        QuESo_ASSERT( p_key->Index() < mSubDictionaries.size(), "It seems like, we forgot to resize 'mSubDictionaries'.\n");

        QuESo_ERROR_IF( !mSubDictionaries[p_key->Index()] )
            << "Dictionary associated with Key: " << std::string( p_key->Name() ) << " is not set.\n";

        return *(mSubDictionaries[p_key->Index()]);
    }

    /// @brief Returns List associated with the given KeyName.
    /// @param rQueryKeyName (std::string).
    /// @return ListType&
    ListType& GetList(const std::string& rQueryKeyName) {

        QuESo_ERROR_IF( !mpListKeySetInfo ) << "This dictionary does not contain any lists.\n";

        const auto p_key = mpListKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName
                               << "'. Possible names are: " << mpListKeySetInfo->GetAllKeyNames() << "\n.";

        QuESo_ASSERT( p_key->Index() < mListsOfDicts.size(), "It seems like, we forgot to resize 'mListsOfDicts'.\n");

        return mListsOfDicts[p_key->Index()];
    }

    /// @brief Sets the value associated with the given KeyName.
    /// @tparam TValueType
    /// @param rQueryKeyName (std::string)
    /// @param rNewValue
    template<typename TValueType>
    void SetValue(const std::string& rQueryKeyName, const TValueType& rNewValue) {
        QuESo_ERROR_IF( !mpDataSet ) << "This dictionary contains an empty data set. You can not set any values here.\n";
        return DataStringAccess::template SetValue<TValueType>(*mpDataSet, rQueryKeyName, rNewValue);
    }

    /// @brief Returns the value associated with the given KeyName.
    /// @tparam TValueType
    /// @param rQueryKeyName (std::string).
    /// @return const TValueType&
    template<typename TValueType>
    const TValueType& GetValue(const std::string& rQueryKeyName) const {
        QuESo_ERROR_IF( !mpDataSet ) << "This dictionary contains an empty data set.\n";
        return DataStringAccess::template GetValue<TValueType>(*mpDataSet, rQueryKeyName);
    }

    /// @brief Returns true if the value associated with the given KeyName is set.
    /// @param rQueryKeyName (std::string).
    /// @return bool
    bool IsSet(const std::string& rQueryKeyName) const {
        QuESo_ERROR_IF( !mpDataSet ) << "This dictionary contains an empty data set.\n";
        return DataStringAccess::IsSet(*mpDataSet, rQueryKeyName);
    }

    /// @brief Returns true if the given KeyName exists.
    /// @param rQueryKeyName (std::string)
    /// @return bool
    bool Has(const std::string& rQueryKeyName) const noexcept {
        const bool has = (mpDataSet && DataStringAccess::Has(*mpDataSet, rQueryKeyName) )
                         || (mpSubDictKeySetInfo && (mpSubDictKeySetInfo->pGetKey(rQueryKeyName) != nullptr) )
                         || (mpListKeySetInfo && (mpListKeySetInfo->pGetKey(rQueryKeyName) != nullptr) );
        return has;
    }

    ///@}
    ///@name Private member variables.
    ///@{

    /// DataSet (Stores all values).
    Unique<DataSet<KeySetValuesTypeTag>> mpDataSet = nullptr;

    /// Members related to SubDictionaries.
    std::vector<Unique<Dictionary>> mSubDictionaries;
    Unique<key::detail::KeySetInfo> mpSubDictKeySetInfo = nullptr;

    /// Members related to Lists.
    std::vector<ListType> mListsOfDicts;
    Unique<key::detail::KeySetInfo> mpListKeySetInfo = nullptr;

    ///@}
}; // End Dictionary class.

/// Output stream function for Dictionary.
template<typename Tag>
inline std::ostream& operator<< (std::ostream& rOStream, const Dictionary<Tag>& rDictionary) {
    std::string indent = "";
    rOStream << std::setprecision(10);
    rOStream << "{\n";
    rDictionary.PrintInfo(rOStream, indent);
    rOStream << '}';
    return rOStream;
}

/// @class DictionaryStringAccess (FOR USE IN PYTHON ONLY).
/// @brief Allows to access the data of Dictionary via KeyNames (std::strings).
///        The respective members are private within Dictionary and made public here.
/// @tparam DictionaryType
/// @details DictionaryStringAccess is a friend of Dictionary
template<typename DictionaryType>
class DictionaryStringAccess {
public:
    ///@name Static operations
    ///@{

    /// @brief Wrapper for Dictionary::SetSubDictionary.
    /// @param rDictionary
    /// @param rQueryKeyName (std::string)
    /// @param rNewDictionary
    static void SetSubDictionary(DictionaryType& rDictionary, const std::string& rQueryKeyName, Unique<DictionaryType>&& rNewDictionary) {
        rDictionary.SetSubDictionary(rQueryKeyName, std::move(rNewDictionary));
    }

    /// @brief Wrapper for Dictionary::SetValue.
    /// @tparam TValueType
    /// @param rDictionary
    /// @param rQueryKeyName (std::string)
    /// @param rNewValue
    template<typename TValueType>
    static void SetValue(DictionaryType& rDictionary, const std::string& rQueryKeyName, const TValueType& rNewValue) {
        rDictionary.SetValue(rQueryKeyName, rNewValue);
    }

    /// @brief Wrapper for Dictionary::GetValue.
    /// @tparam TValueType
    /// @param rDictionary
    /// @param rQueryKeyName (std::string)
    /// @return const TValueType&
    template<typename TValueType>
    static const TValueType& GetValue(const DictionaryType& rDictionary, const std::string& rQueryKeyName) {
        return rDictionary.template GetValue<TValueType>(rQueryKeyName);
    }

    /// @brief Wrapper for Dictionary::GetList.
    /// @param rDictionary
    /// @param rQueryKeyName
    /// @return ListType&
    static typename DictionaryType::ListType& GetList(DictionaryType& rDictionary, const std::string& rQueryKeyName) {
        return rDictionary.GetList(rQueryKeyName);
    }

    /// @brief Wrapper for Dictionary::operator[].
    /// @param rDictionary
    /// @param rQueryKeyName
    /// @return ListType&
    static DictionaryType& GetSubDictionary(DictionaryType& rDictionary, const std::string& rQueryKeyName) {
        return rDictionary[rQueryKeyName];
    }

    /// @brief Wrapper for Dictionary::IsSet.
    /// @param rDictionary
    /// @param rQueryKeyName (std::string)
    /// @return bool
    static bool IsSet(const DictionaryType& rDictionary, const std::string& rQueryKeyName) {
        return rDictionary.IsSet(rQueryKeyName);
    }

    /// @brief Wrapper for Dictionary::Has.
    /// @param rDictionary
    /// @param rQueryKeyName (std::string)
    /// @return bool
    static bool Has(const DictionaryType& rDictionary, const std::string& rQueryKeyName) noexcept {
        return rDictionary.Has(rQueryKeyName);
    }
    ///@}

}; // End class DictionaryStringAccess<>

///@} // End QuESo Classes
} // End queso namespace.

#endif // End DICTIONARY_INCLUDE_HPP
