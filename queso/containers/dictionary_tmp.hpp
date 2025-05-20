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

    // /// @brief Prints this dictionary in JSON format.
    // /// @param rOStream
    // /// @note Precision of doubles is set to 10.
    // void PrintInfo(std::ostream& rOStream) const {
    //     std::string indent = "";
    //     rOStream << std::setprecision(10);
    //     PrintInfo(rOStream, indent);
    // }

private:

    // /// @brief Prints this dictionary in JSON format.
    // /// @param rOStream
    // /// @param rIndent This variable is used for recursive function calls.
    // /// @param IsLastItem This variable is used for recursive function calls.
    // void PrintInfo(std::ostream& rOStream, std::string& rIndent, bool IsLastItem=false) const {
    //     /// Add proper header (with old indent).
    //     if( rIndent == "" ){
    //         rOStream << "{\n"; // Start root
    //     } else if (mKey.index() == 0) { // Key is in std::monostate
    //         rOStream << rIndent << "{\n"; // Start list item
    //     } else if (mIsList) {
    //         rOStream << rIndent << '\"' << mKeyName << "\" : [\n"; // Start List
    //     } else {
    //         rOStream << rIndent << '\"' << mKeyName << "\" : {\n"; // Start dictionary
    //     }

    //     std::string new_indent = rIndent + "\t";
    //     /// Write KeyValuePairs with new indent.
    //     for( IndexType i = 0; i < mData.size(); ++i){
    //         const auto& r_key_value_pair = mData[i];
    //         rOStream << new_indent << "\"" << r_key_value_pair.GetKeyName() << "\" : ";
    //         r_key_value_pair.PrintValue(rOStream);
    //         if( i == mData.size() - 1 && mSubDictionaries.size() == 0 && mListsOfDicts.size() == 0){
    //             rOStream << '\n';
    //         } else {
    //             rOStream << ",\n";
    //         }
    //     }

    //     /// Write mSubDictionaries with new indent.
    //     for( IndexType i = 0; i < mSubDictionaries.size(); ++i){
    //         if( i == mSubDictionaries.size() - 1 && mListsOfDicts.size() == 0) {
    //             mSubDictionaries[i].PrintInfo(rOStream, new_indent, true);
    //         } else {
    //             mSubDictionaries[i].PrintInfo(rOStream, new_indent, false);
    //         }
    //     }

    //     /// Write mListsOfDicts with new indent.
    //     for( IndexType i = 0; i < mListsOfDicts.size(); ++i){
    //         if( i == mListsOfDicts.size() - 1) {
    //             mListsOfDicts[i].PrintInfo(rOStream, new_indent, true);
    //         } else {
    //             mListsOfDicts[i].PrintInfo(rOStream, new_indent, false);
    //         }
    //     }

    //     /// Close with old indent.
    //     if( rIndent == "" || IsLastItem ){ // End list or dict item.
    //         rOStream << ((mIsList) ? (rIndent + "]") : (rIndent + "}")) << "\n"; // Needs to be encapsulated!
    //     } else {
    //         rOStream << ((mIsList) ? (rIndent + "],") : (rIndent + "},")) << "\n";
    //     }
    // }

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
// template<typename... TEnumKeys>
// inline std::ostream& operator<< (std::ostream& rOStream, const Dictionary<TEnumKeys ...>& rThis) {
//     rThis.PrintInfo(rOStream);
//     return rOStream;
// }
///@} // End QuESo Classes
} // End queso namespace.

#endif // End DICTIONARY_INCLUDE_HPP