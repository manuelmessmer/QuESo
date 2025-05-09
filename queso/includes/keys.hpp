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
//

#ifndef KEYS_INCLUDE_HPP
#define KEYS_INCLUDE_HPP

// When Keys are released: remove this include!!
#include "queso/includes/define.hpp"

/// STL includes
#include <string>
#include <vector>
#include <unordered_map>
#include <typeindex>

namespace queso {
namespace key {
namespace detail {

// Count how many commas in the string (i.e., number of items = commas + 1)
constexpr IndexType CountItemsDelimitedByComma(std::string_view str) {
    IndexType count = 1;
    for (char c : str) {
        if (c == ',') ++count;
    }
    return count;
}

// Split a comma-separated string_view into a constexpr std::array
template<IndexType N>
constexpr std::array<std::string_view, N> SplitConstexprStringView(std::string_view str) {
    std::array<std::string_view, N> result{};
    IndexType start = 0;
    IndexType index = 0;

    for (IndexType i = 0; i <= str.size(); ++i) {
        if( str.data()[i] == ' ' || str.data()[i] == '\t' ){ // Skip spaces
            start = i + 1;
        }
        else if (i == str.size() || str.data()[i] == ',') {
            result[index++] = str.substr(start, i - start);
            start = i + 1;
        }
    }

    return result;
}

// Helper to combine counting and splitting
template<IndexType N, const char* TStr>
constexpr auto CreateConstexprStringArray() {
    return SplitConstexprStringView<N>(TStr);
}

/// @brief Helper struct: Allows to define a list of types.
template<typename... Types>
struct ValuesTypeList {
};

/// Forward declaration of KeyToWhatTypeTag template.
template<typename TDerived, typename TTypeList>
struct KeyToWhatTypeTag;

/// @brief This class acts as a type tag to be used for "static_assert". It allows to check to what a
///        specific key leads to.
/// @tparam TDerived - CRTP: Curiously Recurring Template Pattern. This must be the derived class.
/// @tparam TTypeList - List of possible types for the underyling values. Example: You may create a
///         KeyToValues type tag that allows int and double as underlying values: -> ValuesTypeList<int, double>
///         However, you can also create a KeyToSubDict type tag. Since the key always leads to a SubDict,
///         The ValuesTypeList must be kept empty -> ValuesTypeList<>.
template<typename TDerived, typename... PossibleValueTypes>
struct KeyToWhatTypeTag<TDerived, ValuesTypeList<PossibleValueTypes...>> {
    // This is a pure static class.
    // Delete constructor
    KeyToWhatTypeTag() = delete;

    // Delete copy/move constructor and assignment operator
    KeyToWhatTypeTag(const KeyToWhatTypeTag&) = delete;
    KeyToWhatTypeTag& operator=(const KeyToWhatTypeTag&) = delete;
    KeyToWhatTypeTag(KeyToWhatTypeTag&&) = delete;
    KeyToWhatTypeTag& operator=(KeyToWhatTypeTag&&) = delete;

    // Type traits to check, if a given type is a possible value type.
    template<typename TType>
    using is_valid_value_type = std::disjunction<std::is_same<TType, PossibleValueTypes>...>;

    template<typename TType>
    static constexpr bool is_valid_value_type_v = is_valid_value_type<TType>::value;

    /// @brief Returns the name of a given type.
    /// @return std::string
    static constexpr std::string_view GetTypeTagName() {
        return TDerived::msTypeName;
    }

    /// @brief Returns the name of a given type.
    /// @tparam TType
    /// @return std::string
    template<typename TType>
    static constexpr std::string_view GetBaseTypeName() {
        if constexpr(KeyToWhatTypeTag::is_valid_value_type_v<TType>){ // To print name of values.
            return ValuesTypeList<PossibleValueTypes...>::msTypeNames[IndexOfType<TType, PossibleValueTypes...>::msIndex];
        }
        else {
            static_assert(KeyToWhatTypeTag::is_valid_value_type_v<TType>, "Invalid type provided to GetTypeName");
            return "";
        }
    }

private:
    // Primary template, which is not defined.
    template <typename T, typename... Types>
    struct IndexOfType;

    // Specialization 1: T matches the first type in the pack.
    template <typename T, typename... Types>
    struct IndexOfType<T, T, Types...> {
        static constexpr IndexType msIndex = 0;
    };

    // Specialization 2: T does not match the first type in the pack.
    template <typename T, typename U, typename... Types>
    struct IndexOfType<T, U, Types...> {
        static constexpr IndexType msIndex = 1 + IndexOfType<T, Types...>::msIndex;
    };

    // Specialization 3: Parameter pack is empty.
    template <typename T>
    struct IndexOfType<T> {
        static_assert(sizeof(T) == 0, "Type not found in type list!");
    };
};

// Non-template base struct with just the data
template<typename TEnumType>
struct KeyData {
    constexpr KeyData(TEnumType KeyValue) : mKeyValue(KeyValue) {}
    TEnumType mKeyValue;
};

/// Base class of keys.
struct KeyBase {
    virtual ~KeyBase() = default;

    virtual IndexType Index() const = 0;
    virtual std::string_view Name() const = 0;
    virtual std::type_index VariableTypeIndex() const = 0;
    virtual std::type_index KeySetInfoTypeIndex() const = 0;
};

/// @brief Key class to define global keys.
/// @tparam TKeySetInfoType
/// @tparam TKeyToWhat
template<typename TKeySetInfoType, typename TKeyToWhat = typename TKeySetInfoType::KeySetToWhat>
struct KeyConstexpr : private KeyData<typename TKeySetInfoType::EnumType> {
    /// The following typedefs can be used for static type assertion.
    typedef TKeySetInfoType KeySetInfoType;
    typedef TKeyToWhat KeyToWhat;
    typedef typename KeySetInfoType::EnumType KeyValueType;

    /// @brief  Constructor.
    /// @param KeyValue
    constexpr KeyConstexpr(KeyValueType KeyValue) : KeyData<KeyValueType>(KeyValue) {
    }

    /// @brief Returns value/index of Key.
    /// @return IndexType
    constexpr IndexType Index() const noexcept {
        return static_cast<IndexType>(this->mKeyValue);
    }

    /// @brief Returns name of Key.
    /// @return const std::string&
    constexpr std::string_view Name() const noexcept {
        return KeySetInfoType::msEnumNames[static_cast<IndexType>(this->mKeyValue)];
    }
};


/// @brief Key class to define global keys.
/// @tparam TKeySetInfoType
/// @tparam TKeyToWhat
template<typename TKeySetInfoType, typename TKeyToWhat = typename TKeySetInfoType::KeySetToWhat>
struct Key : public KeyBase, private KeyData<typename TKeySetInfoType::EnumType> {
    // Typedefs
    typedef TKeySetInfoType KeySetInfoType;
    typedef TKeyToWhat KeyToWhat;
    typedef typename KeySetInfoType::EnumType KeyValueType;

    /// @brief  Constructor.
    /// @param KeyValue
    Key(KeyValueType KeyValue) : KeyData<KeyValueType>(KeyValue) {
    }

    /// @brief Returns value/index of Key.
    /// @return IndexType
    IndexType Index() const override{
        return static_cast<IndexType>(this->mKeyValue);
    }

    /// @brief Returns name of Key.
    /// @return const std::string&
    std::string_view Name() const override {
        return KeySetInfoType::msEnumNames[static_cast<IndexType>(this->mKeyValue)];
    }

    /// The following functions can be used for dynamic type checks.

    /// @brief Returns std::type_index of type this key points to.
    /// @return std::type_index
    std::type_index VariableTypeIndex() const override {
        return std::type_index(typeid(KeyToWhat));
    }

    /// @brief  Returns std::type_index of type of the key set info, to which this key belongs.
    /// @return std::type_index
    std::type_index KeySetInfoTypeIndex() const override {
        return std::type_index(typeid(KeySetInfoType));
    }

};

/// Base class to store and access information about a key set.
struct KeySetInfo {
    virtual ~KeySetInfo() = default;

    virtual const KeyBase* pGetKey(std::string_view rName) const = 0;
    virtual IndexType GetNumberOfKeys() const noexcept = 0;
    virtual bool IsCorrectKeyType(const KeyBase& rKey) const noexcept = 0;
};

} // End namespace detail
} // End namespace key
} // End namespace queso

#define QuESo_CREATE_VALUE_TYPE_LIST(ValuesTypeListName, ...) \
namespace key {\
    template<> \
    struct detail::ValuesTypeList<__VA_ARGS__> { \
    private:\
        static constexpr char msTypeNamesRaw[] = #__VA_ARGS__; \
        static constexpr std::size_t msSize = queso::key::detail::CountItemsDelimitedByComma(msTypeNamesRaw); \
    public:\
        static constexpr std::array<std::string_view, msSize> msTypeNames = queso::key::detail::CreateConstexprStringArray<msSize, msTypeNamesRaw>();\
    }; \
    typedef detail::ValuesTypeList<__VA_ARGS__> ValuesTypeListName; \
}\

#define QuESo_CREATE_KEY_SET_TO_OBJECT_TYPE_TAG(TypeTagName) \
namespace key {\
    struct TypeTagName : public detail::KeyToWhatTypeTag<TypeTagName, detail::ValuesTypeList<>> {\
        static constexpr char msTypeName[] = #TypeTagName;\
    };\
}\

#define QuESo_CREATE_KEY_SET_TO_VALUE_TYPE_TAG(TypeTagName, ValueTypeList_) \
namespace key {\
    struct TypeTagName : public detail::KeyToWhatTypeTag<TypeTagName, ValueTypeList_> {\
        static constexpr char msTypeName[] = #TypeTagName;\
    };\
}\

#define QuESo_KEY_LIST(...) __VA_ARGS__

#define QuESo_CREATE_KEY_SET_INFO(KeySetName, KeySetToWhat_, ...) \
    typedef queso::key::detail::KeySetInfo KeySetInfoType;\
    typedef queso::key::detail::KeyBase KeyBaseType;\
    namespace key {\
    namespace detail {\
        /* KeySetName##KeyType##KeySetInfo allos to access key set information. */\
        struct KeySetName##KeySetToWhat_##KeySetInfo : public KeySetInfoType {\
            typedef queso::key::KeySetToWhat_ KeySetToWhat;\
            enum class EnumType {__VA_ARGS__};\
            /* Member function to access key information*/\
            const KeyBaseType* pGetKey(std::string_view rName) const override {\
                const auto it = msStringToKeyMap.find(rName);\
                if (it != msStringToKeyMap.end()) {\
                    return (it->second);\
                }\
                QuESo_ERROR << "Invalid Key name. Possible names are: " + StaticGetAllKeyNames();\
            }\
            IndexType GetNumberOfKeys() const noexcept override {\
                return msNumOfEnums;\
            }\
            bool IsCorrectKeyType(const KeyBaseType& rKey) const noexcept override {\
                return (std::type_index(typeid(KeySetName##KeySetToWhat_##KeySetInfo)) == rKey.KeySetInfoTypeIndex() );\
            }\
            static std::string StaticGetAllKeyNames() {\
                std::string result;\
                for (const auto& r_name : msEnumNames) {\
                    if (!result.empty()) {\
                        result += "', '";\
                    }\
                    result += r_name;\
                }\
                return result.insert(0, "['") + "']";\
            }\
        private:\
            static constexpr char msEnumNamesRaw[] = #__VA_ARGS__; \
            static constexpr std::size_t msNumOfEnums = queso::key::detail::CountItemsDelimitedByComma(msEnumNamesRaw); \
        public:\
            static constexpr std::array<std::string_view, msNumOfEnums> msEnumNames = queso::key::detail::CreateConstexprStringArray<msNumOfEnums, msEnumNamesRaw>();\
            static const queso::key::detail::StringToKeyMapType msStringToKeyMap; \
        };\
    }\
    }\

#define QuESo_DEFINE_KEY_SET(KeySetName, KeySetToWhat, KeyNames) \
    QuESo_CREATE_KEY_SET_INFO(KeySetName, KeySetToWhat, KeyNames) \

#define QuESo_DEFINE_KEY_TO_OBJECT(KeySetName, KeyName, KeySetToWhat) \
namespace KeySetName {\
    inline const queso::key::detail::Key<key::detail::KeySetName##KeySetToWhat##KeySetInfo> KeyName(\
        key::detail::KeySetName##KeySetToWhat##KeySetInfo::EnumType::KeyName);\
}\

#define QuESo_DEFINE_KEY_TO_VALUE(KeySetName, KeyName, KeySetToWhat, KeyToWhat) \
    namespace KeySetName {\
        static_assert(queso::key::KeySetToWhat::is_valid_value_type_v<KeyToWhat>, "Given KeyToWhat-type is invalid.");\
        inline const queso::key::detail::Key<key::detail::KeySetName##KeySetToWhat##KeySetInfo, KeyToWhat> KeyName(\
            key::detail::KeySetName##KeySetToWhat##KeySetInfo::EnumType::KeyName);\
        inline constexpr queso::key::detail::KeyConstexpr<key::detail::KeySetName##KeySetToWhat##KeySetInfo, KeyToWhat> KeyName##2(\
            key::detail::KeySetName##KeySetToWhat##KeySetInfo::EnumType::KeyName);\
    }\

namespace queso {
namespace key {
namespace detail {

typedef std::unordered_map<std::string_view, const queso::key::detail::KeyBase* const> StringToKeyMapType;

/// @brief Helper functions to check if all keys are correctly defined.
/// @tparam TKeySetInfoType
/// @param rKeyMap
/// @return StringToKeyMapType
/// @note This function is called during the runtime initialization of the respective static object.
template<typename TKeySetInfoType>
inline StringToKeyMapType InitializeStringToKeyMap(const StringToKeyMapType& rKeyMap) {
    if( rKeyMap.size() != TKeySetInfoType::msEnumNames.size() ) {
        QuESo_ERROR << "Number of declared keys does not match number of defined keys. Check: REGISTER_KEYS_HPP.";
    }
    for( const auto& r_key_name : TKeySetInfoType::msEnumNames ) {
        const auto it = rKeyMap.find(r_key_name);
        if (it == rKeyMap.end()) {
            QuESo_ERROR << "Key '" << r_key_name << "' not found in key map. Check: REGISTER_KEYS_HPP.";
        }
    }
    return rKeyMap;
}

} // End namespace detail
} // End namespace key
} // End namespace queso

#define QuESo_KEY(Key) {Key.Name(), &Key} \

#define QuESo_REGISTER_KEY_SET(KeySet, KeySetToWhat, ...) \
    inline const queso::key::detail::StringToKeyMapType key::detail::KeySet##KeySetToWhat##KeySetInfo::msStringToKeyMap =\
        queso::key::detail::InitializeStringToKeyMap<key::detail::KeySet##KeySetToWhat##KeySetInfo>( {__VA_ARGS__} );\


#endif // End KEYS_INCLUDE_HPP