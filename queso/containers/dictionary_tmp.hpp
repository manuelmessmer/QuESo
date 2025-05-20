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
#include "queso/containers/data_set.hpp"
#include "queso/includes/define.hpp"

namespace queso {

/// Forward declaration
template<typename DictionaryType>
class DictionaryStringAccess;

///@name QuESo Classes
///@{

/**
 * @class  Dictionary.
 * @author Manuel Messmer
 * @brief Base class to store json-like data sets. We can store KeyValuePair's, Subdictionaries, or Lists.
 *        The Dictionary is optimized to allow extremely fast access to the underlying data. Therefore, the keys are hardcoded.
 *        The KeyTypes are supposed to be enum classes and can be passed as a template parameter pack: ...TEnumKeys.
 *        The main idea is to derive from Dictionary and to 'hardcode' the respective keys/enums.
 *        Thereby, the keys on the same level, must be of the same type (same enum). ( @see e.g. Settings ).
 *        Most functions also allow to pass the Key as and std::string. However, this is intended only to be used in Python.
**/
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

    // Helper type tag to use a constructor template.
    template<typename TValueKeySetInfoType = std::false_type,
             typename TSubDictKeySetInfoType = std::false_type,
             typename TListKeySetInfoType = std::false_type>
    struct KeySetInfoTypeTag {
        using ValueKeySetInfoType = TValueKeySetInfoType;
        using SubDictKeySetInfoType = TSubDictKeySetInfoType;
        using ListKeySetInfoType = TListKeySetInfoType;
    };

    // Type trait helpers
    template<typename TType>
    static constexpr bool has_value_key_v = !std::is_same_v<typename TType::ValueKeySetInfoType, std::false_type>;

    template<typename TType>
    static constexpr bool has_subdict_key_v = !std::is_same_v<typename TType::SubDictKeySetInfoType, std::false_type>;

    template<typename TType>
    static constexpr bool has_list_key_v = !std::is_same_v<typename TType::ListKeySetInfoType, std::false_type>;

    ///@}
    ///@name Life cycle
    ///@{
    Dictionary() = default;

    template<typename TKeySetInfoTypeTag>
    Dictionary( TKeySetInfoTypeTag ) {
        using Tag = TKeySetInfoTypeTag;

        if constexpr (has_value_key_v<Tag>) {
            using KeySetInfoType = typename Tag::ValueKeySetInfoType;
            using DataSetTypeTag = typename DataSetType::template KeySetInfoTypeTag<KeySetInfoType>;
            mpDataSet = MakeUnique<DataSetType>( DataSetTypeTag{} );
        }
        if constexpr ( has_subdict_key_v<Tag> ) {
            using KeySetInfoType = typename Tag::SubDictKeySetInfoType;
            mpSubDictKeySetInfo = MakeUnique<KeySetInfoType>();
            mSubDictionaries.resize(mpSubDictKeySetInfo->GetNumberOfKeys());
        }
        if constexpr ( has_list_key_v<Tag> ) {
            using KeySetInfoType = typename Tag::ListKeySetInfoType;
            mpListKeySetInfo = MakeUnique<KeySetInfoType>();
            mListsOfDicts.resize(mpListKeySetInfo->GetNumberOfKeys());
        }
    }

    /// Destructor
    virtual ~Dictionary() = default;
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

    template<typename TKeyType, typename TValueType>
    void SetValue(const TKeyType& rQueryKey, const TValueType& rNewValue) noexcept(NOTDEBUG) {
        QuESo_ASSERT( mpDataSet != nullptr, "This dictionary has an empty data set.\n");
        mpDataSet->SetValue(rQueryKey, rNewValue);
    }

    template<typename TValueType, typename TKeyType>
    const TValueType& GetValueFast(const TKeyType& rQueryKey) const noexcept(NOTDEBUG) {
        QuESo_ASSERT( mpDataSet != nullptr, "This dictionary has an empty data set.\n");
        return mpDataSet->template GetValueFast<TValueType>(rQueryKey);
    }

    template<typename TValueType, typename TKeyType>
    const TValueType& GetValue(const TKeyType& rQueryKey) const {
        QuESo_ASSERT( mpDataSet != nullptr, "This dictionary has an empty data set.\n");
        return mpDataSet->template GetValue<TValueType>(rQueryKey);
    }

    template<typename TKeyType>
    bool IsSet(const TKeyType& rQueryKey) const noexcept(NOTDEBUG) {
        QuESo_ASSERT( mpDataSet != nullptr, "This dictionary has an empty data set.\n");
        return mpDataSet->IsSet(rQueryKey);
    }

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

    /// @brief Returns ref to subdictionary with the given Key (const version).
    /// @tparam TKeyType.s
    /// @param rQueryKey
    /// @return const Dictionary&
    /// @note only throws in DEBUG mode.
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

    /// @brief Returns ref to subdictionary corresponding to given key provided as Enum (non const version).
    /// @tparam TKeyType.
    /// @param rQueryKey
    /// @return Dictionary&
    /// @note only throws in DEBUG mode.
    template<typename TKeyType,
             typename = std::enable_if_t<is_key_v<TKeyType>>>
    Dictionary& operator[](const TKeyType& rQueryKey) noexcept(NOTDEBUG) {
        return const_cast<Dictionary&>( static_cast<const Dictionary&>(*this)[rQueryKey] );
    }

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

    template<typename TKeyType>
    ListType& GetList(const TKeyType& rQueryKey,
                      std::enable_if_t<is_key_v<TKeyType>>* = nullptr  ) noexcept(NOTDEBUG) {
        return const_cast<ListType&>(static_cast<const Dictionary&>(*this).GetList(rQueryKey));
    }

    /// @brief Prints this dictionary in JSON format.
    /// @param rOStream
    /// @note Precision of doubles is set to 10.
    void PrintInfo(std::ostream& rOStream) const {
        std::string indent = "";
        rOStream << std::setprecision(10);
        rOStream << "{\n";
        PrintInfo(rOStream, indent);
        rOStream << '}';
    }

private:
    ///@}
    ///@name Type definitions
    ///@{

    using DataStringAccess = DataSetStringAccess<DataSetType>;

    ///@}
    ///@name Private operations
    ///@{

    /// @brief Prints this dictionary in JSON format.
    /// @param rOStream
    /// @param rIndent This variable is used for recursive function calls.
    /// @param IsLastItem This variable is used for recursive function calls.
    void PrintInfo(std::ostream& rOStream, std::string& rIndent, bool IsLastItem=false) const {
        std::string new_indent = rIndent + "\t";
        if( mpDataSet ) {
            const auto& r_data_key_set_info = mpDataSet->GetKeySetInfo();
            const auto& r_map = r_data_key_set_info.GetStringKeyMap();
            auto it = r_map.begin();

            while (it != r_map.end()) {
                const auto& [name, p_key] = *it;

                rOStream << new_indent << "\"" << name << "\" : ";
                mpDataSet->PrintValue(p_key, rOStream);

                if (++it == r_map.end()) {
                    rOStream << '\n'; // last element
                } else {
                    rOStream << ",\n";
                }
            }
        }

        /// Write mSubDictionaries with new indent.
        if( mpSubDictKeySetInfo ) {
            const auto& r_map = mpSubDictKeySetInfo->GetStringKeyMap();
            auto it = r_map.begin();

            while (it != r_map.end()) {
                const auto& [r_name, p_key] = *it;
                const auto& p_sub_dict = mSubDictionaries[p_key->Index()];
                ++it;
                if( p_sub_dict ) {
                    rOStream << new_indent << '\"' << r_name << "\" : {\n"; // Start dictionary
                    p_sub_dict->PrintInfo(rOStream, new_indent, true);
                    if( it == r_map.end() && mListsOfDicts.size() == 0) {
                        rOStream << (new_indent + "}\n");
                    } else {
                        rOStream << (new_indent + "},\n");
                    }
                }
            }
        }

        /// Write mSubDictionaries with new indent.
        if( mpListKeySetInfo ) {
            const auto& r_map = mpListKeySetInfo->GetStringKeyMap();
            auto it = r_map.begin();

            while (it != r_map.end()) {
                const auto& [r_name, p_key] = *it;
                const auto& r_list = mListsOfDicts[p_key->Index()];

                rOStream << new_indent << '\"' << r_name << "\" : [\n"; // Start List
                std::string inner_indent = new_indent + "\t";
                for( IndexType i = 0; i < r_list.size() ; ++i){
                    rOStream << inner_indent << "{\n"; // Start list item
                    r_list[i]->PrintInfo(rOStream, inner_indent, true);
                    if( i == r_list.size() - 1) {
                        rOStream << (inner_indent + "}\n");
                    } else {
                        rOStream << (inner_indent + "},\n");
                    }
                }

                if( ++it == r_map.end() ) {
                    rOStream << (new_indent + "}\n");
                } else {
                    rOStream << (new_indent + "},\n");
                }

            }
        }
    }

    /// The following operations allow to access the data via KeyNames (std::string).
    /// These functions should only be used in Python. Therefore, they are made private here
    /// and can only be used via the friend class DictionaryStringAccess.
    friend class DictionaryStringAccess<Dictionary>;

    void SetSubDictionary(const std::string& rQueryKeyName, Unique<Dictionary>&& rNewDictionary) {

        QuESo_ERROR_IF( !mpSubDictKeySetInfo )  << "This dictionary does not contain any sub-dictionaries.\n";

        const auto p_key = mpSubDictKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName
                               << "'. Possible names are: " << mpSubDictKeySetInfo->GetAllKeyNames() << "\n.";

        QuESo_ASSERT( p_key->Index() < mSubDictionaries.size(), "It seems like, we forgot to resize 'mSubDictionaries'.\n");

        mSubDictionaries[p_key->Index()] = std::move(rNewDictionary);
    }


    /// @brief Returns ref to subdictionary corresponding to given key provided as std::string.
    ///        This function should only be used in Python. The Version: [Enum Key] is much faster.
    /// @param rQueryKeyName
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


    ListType& GetList(const std::string& rQueryKeyName) {

        QuESo_ERROR_IF( !mpListKeySetInfo ) << "This dictionary does not contain any lists.\n";

        const auto p_key = mpListKeySetInfo->pGetKey(rQueryKeyName);
        QuESo_ERROR_IF(!p_key) << "Invalid Key name: '" << rQueryKeyName
                               << "'. Possible names are: " << mpListKeySetInfo->GetAllKeyNames() << "\n.";

        QuESo_ASSERT( p_key->Index() < mListsOfDicts.size(), "It seems like, we forgot to resize 'mListsOfDicts'.\n");

        return mListsOfDicts[p_key->Index()];
    }

    template<typename TValueType>
    void SetValue(const std::string& rQueryKeyName, const TValueType& rNewValue) {
        QuESo_ERROR_IF( !mpDataSet ) << "This dictionary contains an empty data set. You can not set any values here.\n";
        return DataStringAccess::template SetValue<TValueType>(*mpDataSet, rQueryKeyName, rNewValue);
    }

    template<typename TValueType>
    const TValueType& GetValue(const std::string& rQueryKeyName) const {
        QuESo_ERROR_IF( !mpDataSet ) << "This dictionary contains an empty data set.\n";
        return DataStringAccess::template GetValue<TValueType>(*mpDataSet, rQueryKeyName);
    }

    bool IsSet(const std::string& rQueryKeyName) const {
        QuESo_ERROR_IF( !mpDataSet ) << "This dictionary contains an empty data set.\n";
        return DataStringAccess::IsSet(*mpDataSet, rQueryKeyName);
    }

    bool Has(const std::string& rQueryKeyName) const noexcept {
        const bool has = (mpDataSet && DataStringAccess::Has(*mpDataSet, rQueryKeyName) )
                         || (mpSubDictKeySetInfo && (mpSubDictKeySetInfo->pGetKey(rQueryKeyName) != nullptr) )
                         || (mpListKeySetInfo && (mpListKeySetInfo->pGetKey(rQueryKeyName) != nullptr) );
        return has;
    }



    ///@}
    ///@name Private member variables.
    ///@{


    /// Vector of KeyValuePairs.
    Unique<DataSet<KeySetValuesTypeTag>> mpDataSet = nullptr;

    /// Vector of Subdictionaries.
    std::vector<Unique<Dictionary>> mSubDictionaries;
    Unique<key::detail::KeySetInfo> mpSubDictKeySetInfo = nullptr;

    /// Members related to Lists.
    std::vector<ListType> mListsOfDicts;
    Unique<key::detail::KeySetInfo> mpListKeySetInfo = nullptr;

    ///@}
};

/// Output stream functions
/// Make friend
template<typename... TEnumKeys>
inline std::ostream& operator<< (std::ostream& rOStream, const Dictionary<TEnumKeys ...>& rThis) {
    rThis.PrintInfo(rOStream);
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

    /// @brief Wrapper for Dictionary::SetSubDictionary
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

    /// @brief Wrapper for Dictionary[].
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
    /// @param rDictionary rDictionary
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