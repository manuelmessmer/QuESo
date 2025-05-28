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

/// STL includes
#include <variant>
#include <string>
#include <unordered_map>
#include <typeindex>

namespace queso {
namespace key {
namespace detail {

/// @brief Returns the number of items in a string delimited by a comma
///        (i.e., number of items = commas + 1).
/// @param StrView
/// @return IndexType
constexpr IndexType NumberOfItemsDelimitedByComma(std::string_view StrView) {
    IndexType count = 1;
    for (char c : StrView) {
        if (c == ',') ++count;
    }
    return count;
}

/// @brief Splits a string into an array of strings using ',' as delimiter.
/// @tparam TSize
/// @param StrView
/// @return std::array<std::string_view, TSize>
template<IndexType TSize>
constexpr std::array<std::string_view, TSize> CreateConstexprStringArray(std::string_view StrView) {
    std::array<std::string_view, TSize> result{};
    IndexType start = 0;
    IndexType index = 0;

    for (IndexType i = 0; i <= StrView.size(); ++i) {
        if( StrView.data()[i] == ' ' || StrView.data()[i] == '\t' ){ // Skip spaces and tabs.
            start = i + 1;
        }
        else if (i == StrView.size() || StrView.data()[i] == ',') { // Add substring.
            result[index++] = StrView.substr(start, i - start);
            start = i + 1;
        }
    }
    return result;
}

/// @brief Helper struct: Allows to define a list of types.
template<typename... Types>
struct ValuesTypeList {
};

/// Forward declaration of KeyToWhatTypeTag template.
template<typename TDerived, typename TTypeList>
struct KeyToWhatTypeTag;

/// @brief This class acts as a type tag to be used for "static_assert". It allows to check to what type a
///        specific key leads to.
/// @tparam TDerived - CRTP: Curiously Recurring Template Pattern. This must be the derived class.
///         TDerived is the actual type tag, e.g., KeyToExampleValuesTypeTag or KeyToSubDictTypeTag.
/// @tparam TTypeList - List of possible types for the underyling values (may be empty). Example: You may create a
///         KeyToExampleValuesTypeTag that allows 'int' and 'double' as underlying values:
///             - TDerived: KeyToExampleValuesTypeTag
///             - TTypeList: ValuesTypeList<int, double>
///         However, you can also create a KeyToSubDict type tag, where keys always lead/point to another subdict.
///             - TDerived: KeyToSubDictTypeTag
///             - TTypeList: ValuesTypeList<>
template<typename TDerived, typename... PossibleValueTypes>
struct KeyToWhatTypeTag<TDerived, ValuesTypeList<PossibleValueTypes...>> {
    /// This is a pure static class.
    /// Delete constructor.
    KeyToWhatTypeTag() = delete;

    /// Delete copy/move constructor and assignment operators.
    KeyToWhatTypeTag(const KeyToWhatTypeTag&) = delete;
    KeyToWhatTypeTag& operator=(const KeyToWhatTypeTag&) = delete;
    KeyToWhatTypeTag(KeyToWhatTypeTag&&) = delete;
    KeyToWhatTypeTag& operator=(KeyToWhatTypeTag&&) = delete;

    /// Typedef
    using VariantType = std::variant<std::monostate, PossibleValueTypes...>;

    /// Type traits to check, if a given type is a actually a key.
    template<typename, typename = void>
    struct is_key : std::false_type {};

    template<typename TKeyType>
    struct is_key<TKeyType, std::void_t<typename TKeyType::KeyValueType>>
        : std::integral_constant<bool,
            !std::is_convertible_v<typename TKeyType::KeyValueType, int> &&
            std::is_enum_v<typename TKeyType::KeyValueType>> {};

    template<typename TKeyType>
    static constexpr bool is_key_v = is_key<TKeyType>::value;

    /// Type traits to check, if a given type is a possible value type.
    template<typename TType>
    using is_valid_value_type = std::disjunction<std::is_same<TType, PossibleValueTypes>...>;

    template<typename TType>
    static constexpr bool is_valid_value_type_v = is_valid_value_type<TType>::value;

    /// @brief Returns the name of this type tag.
    ///        msTypeName must be specified in derived class.
    /// @return std::string_view
    static constexpr std::string_view GetTypeTagName() noexcept {
        return TDerived::msTypeName;
    }

    /// @brief Returns the name the underlying type (if applicable).
    ///        msTypeNames must be defined in ValuesTypeList<>.
    /// @tparam TType
    /// @return std::string
    template<typename TType>
    static constexpr std::string_view GetValueTypeName() noexcept {
        static_assert(KeyToWhatTypeTag::is_valid_value_type_v<TType>, "Invalid type provided to GetValueTypeName");
        return ValuesTypeList<PossibleValueTypes...>::msTypeNames[IndexOfType<TType, PossibleValueTypes...>::msIndex];
    }

private:
    /// Primary template, which is not defined.
    template <typename T, typename... Types>
    struct IndexOfType;

    /// Specialization 1: T matches the first type in the pack.
    template <typename T, typename... Types>
    struct IndexOfType<T, T, Types...> {
        static constexpr IndexType msIndex = 0;
    };

    /// Specialization 2: T does not match the first type in the pack.
    template <typename T, typename U, typename... Types>
    struct IndexOfType<T, U, Types...> {
        static constexpr IndexType msIndex = 1 + IndexOfType<T, Types...>::msIndex;
    };

    /// Specialization 3: Parameter pack is empty.
    template <typename T>
    struct IndexOfType<T> {
        static_assert(sizeof(T) == 0, "Type not found in type list!");
    };
};

/// @brief  KeyData: stores the KeyValue for all keys.
/// @tparam TEnumType
template<typename TEnumType>
struct KeyData {
    constexpr KeyData(TEnumType KeyValue) : mKeyValue(KeyValue) {}
protected:
    const TEnumType mKeyValue;
};

/// @brief DynamicKeyBase: Base class for dynamic keys.
struct DynamicKeyBase {
    virtual ~DynamicKeyBase() = default;

    /// Pure virtual functions.
    virtual IndexType Index() const noexcept = 0;
    virtual std::string_view Name() const noexcept = 0;
    virtual std::string_view TargetTypeName() const noexcept = 0;
    virtual std::type_index TargetTypeIndex() const noexcept = 0;
    virtual std::type_index KeySetInfoTypeIndex() const noexcept = 0;
};

/// @brief Dynamic key class.
///        Derives publicly from DynamicKeyBase to allow dynamic polymorphism.
///        Derives privately from KeyData to access the respective key value.
/// @tparam TKeySetInfoType
/// @tparam TKeyToWhat (Default: TKeyToWhat = TKeySetInfoType::KeySetToWhat )
template<typename TKeySetInfoType, typename TKeyToWhat = typename TKeySetInfoType::KeySetToWhat>
struct DynamicKey : public DynamicKeyBase, private KeyData<typename TKeySetInfoType::EnumType> {
    /// Typedefs
    using KeySetInfoType = TKeySetInfoType;
    using KeyToWhat = TKeyToWhat;
    using KeyValueType = typename KeySetInfoType::EnumType;

    /// @brief  Constructor.
    /// @param KeyValue
    DynamicKey(KeyValueType KeyValue) : KeyData<KeyValueType>(KeyValue) {
    }

    /// @brief Returns the value/index of this Key.
    /// @return IndexType
    IndexType Index() const noexcept override {
        return static_cast<IndexType>(this->mKeyValue);
    }

    /// @brief Returns the name of this Key.
    /// @return const std::string&
    std::string_view Name() const noexcept override {
        return KeySetInfoType::msEnumNames[static_cast<IndexType>(this->mKeyValue)];
    }

    /// @brief Returns the name of the target type.
    ///        This can be a value, e.g., 'double', 'IndexType', etc. or a KeySetTypeTag, e.g., 'SubDictTypeTag'.
    /// @return std::string_view
    std::string_view TargetTypeName() const noexcept override {
        using KeySetTypeTag = typename KeySetInfoType::KeySetToWhat;
        if constexpr ( KeySetTypeTag::template is_valid_value_type_v<KeyToWhat> ) {
            return KeySetTypeTag::template GetValueTypeName<KeyToWhat>();
        } else {
            return KeySetTypeTag::GetTypeTagName();
        }
    }

    /* The following functions can be used for dynamic type checks. */

    /// @brief Returns std::type_index of type this key points to (KeyToWhat).
    /// @return std::type_index
    std::type_index TargetTypeIndex() const noexcept override {
        return std::type_index(typeid(KeyToWhat));
    }

    /// @brief  Returns std::type_index of type of the KeySetInfo, to which this key belongs.
    /// @return std::type_index
    std::type_index KeySetInfoTypeIndex() const noexcept override {
        return std::type_index(typeid(KeySetInfoType));
    }
};

/// @brief ConstexprKey. All members are constexpr.
///        Derives privately from KeyData to access the respective key value.
/// @tparam TKeySetInfoType
/// @tparam TKeyToWhat
template<typename TKeySetInfoType, typename TKeyToWhat = typename TKeySetInfoType::KeySetToWhat>
struct ConstexprKey : private KeyData<typename TKeySetInfoType::EnumType> {
    /// The following typedefs can be used for static type assertion.
    using KeySetInfoType = TKeySetInfoType;
    using KeyToWhat = TKeyToWhat;
    using KeyValueType = typename KeySetInfoType::EnumType;

    /// @brief  Constructor (constexpr).
    /// @param KeyValue
    constexpr ConstexprKey(KeyValueType KeyValue) : KeyData<KeyValueType>(KeyValue) {
    }

    /// @brief Returns the value/index of this Key.
    /// @return IndexType
    constexpr IndexType Index() const noexcept {
        return static_cast<IndexType>(this->mKeyValue);
    }

    /// @brief Returns the name of this Key.
    /// @return const std::string&
    constexpr std::string_view Name() const noexcept {
        return KeySetInfoType::msEnumNames[static_cast<IndexType>(this->mKeyValue)];
    }

    /// @brief  Returns std::type_index of type of the KeySetInfo, to which this key belongs.
    /// @return std::type_index
    constexpr std::type_index KeySetInfoTypeIndex() const noexcept {
        return std::type_index(typeid(KeySetInfoType));
    }
};

/// Typedef
using StringToKeyMapType = std::unordered_map<std::string_view, const queso::key::detail::DynamicKeyBase*>;

/// Base class to store and access information about a key set.
struct KeySetInfo {
    virtual ~KeySetInfo() = default;

    virtual const DynamicKeyBase* pGetKey(std::string_view rName) const noexcept = 0;
    virtual const StringToKeyMapType& GetStringKeyMap() const noexcept = 0;
    virtual IndexType GetNumberOfKeys() const noexcept = 0;
    virtual bool IsSameKeySet(std::type_index TypeIndex) const noexcept = 0;
    virtual bool IsPartOfKeySet(const DynamicKeyBase& rKey) const noexcept = 0;
    virtual std::string GetAllKeyNames() const = 0;
};

} // End namespace detail
} // End namespace key
} // End namespace queso

/// @brief Macro to define a value type list with associated string names.
/// @param ValuesTypeListName The typedef name for the list.
/// @param ... Variadic list of types (e.g., int, double, std::string).
/// @details Generates a specialization of `ValuesTypeList` that includes
///          a constexpr array of type names as string_view entries for type
///          introspection.
#define QuESo_CREATE_VALUE_TYPE_LIST(ValuesTypeListName, ...)\
namespace key {\
    template<>\
    struct detail::ValuesTypeList<__VA_ARGS__> {\
    private:\
        static constexpr char msTypeNamesRaw[] = #__VA_ARGS__;\
        static constexpr std::size_t msSize = queso::key::detail::NumberOfItemsDelimitedByComma(msTypeNamesRaw);\
    public:\
        static constexpr std::array<std::string_view, msSize> msTypeNames = queso::key::detail::CreateConstexprStringArray<msSize>(msTypeNamesRaw);\
    };\
    using ValuesTypeListName = detail::ValuesTypeList<__VA_ARGS__>;\
}\

/// @brief Macro to define a type tag for keys that point to objects (e.g., dicstionaries or lists).
/// @param TypeTagName Name of the type tag struct.
/// @details Defines a type tag inheriting from `KeyToWhatTypeTag` with an empty value type list.
#define QuESo_CREATE_KEY_SET_TO_OBJECT_TYPE_TAG(TypeTagName)\
namespace key {\
    struct TypeTagName : public detail::KeyToWhatTypeTag<TypeTagName, detail::ValuesTypeList<>> {\
        static constexpr char msTypeName[] = #TypeTagName;\
    };\
}\

/// @brief Macro to define a type tag for keys that point to typed values.
/// @param TypeTagName Name of the type tag struct.
/// @param ValueTypeList The value type list (e.g., created with `QuESo_CREATE_VALUE_TYPE_LIST`).
/// @details Defines a type tag using the specified list of acceptable types.
#define QuESo_CREATE_KEY_SET_TO_VALUE_TYPE_TAG(TypeTagName, ValueTypeList)\
namespace key {\
    struct TypeTagName : public detail::KeyToWhatTypeTag<TypeTagName, ValueTypeList> {\
        static constexpr char msTypeName[] = #TypeTagName;\
    };\
}\

/// @brief Utility macro to simplify writing variadic key lists.
/// @param ... List of key names.
/// @details Meant for syntactic convenience when passing multiple key names.
#define QuESo_KEY_LIST(...) __VA_ARGS__

/// @brief Declares metadata and interfaces for a key set.
/// @param KeySetName Name of the key set.
/// @param KeySetToWhat_ Associated type tag.
/// @param ... List of key names (enum values).
/// @details Defines an enum of keys and a class providing runtime access to key names and validation.
#define QuESo_DECLARE_KEY_SET_INFO(KeySetName, KeySetToWhat_, ...)\
    using KeySetInfoType = queso::key::detail::KeySetInfo;\
    using KeyBaseType = queso::key::detail::DynamicKeyBase;\
    namespace key {\
    namespace detail {\
        struct KeySetName##KeySetToWhat_##KeySetInfo : public KeySetInfoType {\
            using KeySetToWhat = queso::key::KeySetToWhat_;\
            enum class EnumType {__VA_ARGS__};\
            /* Member function to access key information*/\
            const KeyBaseType* pGetKey(std::string_view rName) const noexcept override {\
                const auto it = msStringToKeyMap.find(rName);\
                if (it != msStringToKeyMap.end()) {\
                    return (it->second);\
                }\
                return nullptr;\
            }\
            const queso::key::detail::StringToKeyMapType& GetStringKeyMap() const noexcept override {\
                return msStringToKeyMap;\
            }\
            IndexType GetNumberOfKeys() const noexcept override {\
                return msNumOfEnums;\
            }\
            bool IsSameKeySet(std::type_index TypeIndex) const noexcept override {\
                return (std::type_index(typeid(KeySetName##KeySetToWhat_##KeySetInfo)) == TypeIndex);\
            }\
            bool IsPartOfKeySet(const KeyBaseType& rKey) const noexcept override {\
                return (std::type_index(typeid(KeySetName##KeySetToWhat_##KeySetInfo)) == rKey.KeySetInfoTypeIndex() );\
            }\
            std::string GetAllKeyNames() const override {return StaticGetAllKeyNames(); }\
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
            static constexpr std::size_t msNumOfEnums = queso::key::detail::NumberOfItemsDelimitedByComma(msEnumNamesRaw);\
        public:\
            static constexpr std::array<std::string_view, msNumOfEnums> msEnumNames = queso::key::detail::CreateConstexprStringArray<msNumOfEnums>(msEnumNamesRaw);\
            static const queso::key::detail::StringToKeyMapType msStringToKeyMap;\
        };\
    }\
    }\

/// @brief Convenience macro to define a complete key set with enum and associated metadata.
/// @param KeySetName Name of the key set.
/// @param KeySetToWhat The type tag.
/// @param KeyNames List of keys, passed via `QuESo_KEY_LIST(...)`.
/// @details Internally calls `QuESo_DECLARE_KEY_SET_INFO`.
#define QuESo_DEFINE_KEY_SET(KeySetName, KeySetToWhat, KeyNames)\
    QuESo_DECLARE_KEY_SET_INFO(KeySetName, KeySetToWhat, KeyNames)\

/// @brief Defines a key that points to an object (e.g., dictionary).
/// @param KeySetName Name of the key set.
/// @param KeyName The key identifier (enum entry).
/// @param KeySetToWhat Type tag this key set maps to.
/// @details Creates both a dynamic and a constexpr version of the key.
#define QuESo_DEFINE_KEY_TO_OBJECT(KeySetName, KeyName, KeySetToWhat)\
    namespace key {\
    namespace detail {\
    namespace KeySetName {/* Define dynamic key */\
        inline const queso::key::detail::DynamicKey<key::detail::KeySetName##KeySetToWhat##KeySetInfo> KeyName##_Dynamic(\
            key::detail::KeySetName##KeySetToWhat##KeySetInfo::EnumType::KeyName);\
    }}}\
    namespace KeySetName {/* Define constexpr key */\
        inline constexpr queso::key::detail::ConstexprKey<key::detail::KeySetName##KeySetToWhat##KeySetInfo> KeyName(\
            key::detail::KeySetName##KeySetToWhat##KeySetInfo::EnumType::KeyName);\
    }\

/// @brief Defines a key that points to a value of a specific type.
/// @param KeySetName Name of the key set.
/// @param KeyName The key identifier (enum entry).
/// @param KeySetToWhat The key set's type tag (value category).
/// @param KeyToWhat The specific value type the key maps to (e.g., int, double).
/// @details Generates dynamic and constexpr keys. Includes a static assert to ensure the type is valid.
#define QuESo_DEFINE_KEY_TO_VALUE(KeySetName, KeyName, KeySetToWhat, KeyToWhat)\
    static_assert(queso::key::KeySetToWhat::is_valid_value_type_v<KeyToWhat>, "Given KeyToWhat-type is invalid.");\
    namespace key {\
    namespace detail {\
    namespace KeySetName {/* Define dynamic key */\
        inline const queso::key::detail::DynamicKey<key::detail::KeySetName##KeySetToWhat##KeySetInfo, KeyToWhat> KeyName##_Dynamic(\
            key::detail::KeySetName##KeySetToWhat##KeySetInfo::EnumType::KeyName);\
    }}}\
    namespace KeySetName {/* Define constexpr key */\
        inline constexpr queso::key::detail::ConstexprKey<key::detail::KeySetName##KeySetToWhat##KeySetInfo, KeyToWhat> KeyName(\
            key::detail::KeySetName##KeySetToWhat##KeySetInfo::EnumType::KeyName);\
    }\

namespace queso {
namespace key {
namespace detail {

/// @brief Helper function to check if all keys are correctly defined.
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

/// @brief Defines a key entry for use in the key map registration.
/// @param Key The key name.
/// @details Creates a mapping from string name to pointer to dynamic key object.
#define QuESo_KEY(Key) {key::detail::Key##_Dynamic.Name(), &key::detail::Key##_Dynamic}\

/// @brief Registers a key set's dynamic key map.
/// @param KeySet Name of the key set.
/// @param KeySetToWhat The associated type tag.
/// @param ... List of key mappings created with `QuESo_KEY(...)`.
/// @details Initializes the key map and performs runtime consistency checks.
#define QuESo_REGISTER_KEY_SET(KeySet, KeySetToWhat, ...)\
    inline const queso::key::detail::StringToKeyMapType key::detail::KeySet##KeySetToWhat##KeySetInfo::msStringToKeyMap =\
        queso::key::detail::InitializeStringToKeyMap<key::detail::KeySet##KeySetToWhat##KeySetInfo>( {__VA_ARGS__} );\

#endif // End KEYS_INCLUDE_HPP