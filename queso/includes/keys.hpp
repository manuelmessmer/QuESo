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

#include "queso/includes/define.hpp"
/// STL includes
#include <string>
#include <vector>
#include <unordered_map>
#include <typeindex>

namespace queso {
namespace key {
namespace detail {

typedef std::vector<std::string> StringVectorType;

/// @brief Helper function to remove white spaces from a string.
/// @param rString
/// @return std::string
inline std::string GetStringWithoutWhiteSpaces(const std::string& rString) {
    std::string result;
    result.reserve(rString.size());
    for (char ch : rString) {
        if (!std::isspace(static_cast<unsigned char>(ch))) {
            result.push_back(ch);
        }
    }
    return result;
}

/// @brief Helper function to remove 'QuESo_KEY_LIST(' and ')' from a string.
/// @param rString
inline void RemoveQuESoList(std::string& rString) {
    std::size_t start_pos = rString.find("QuESo_KEY_LIST(");
    if (start_pos != std::string::npos) {
        rString.erase(start_pos, std::string("QuESo_KEY_LIST(").length());
    }
    std::size_t end_pos = rString.find(")", start_pos);
    if (end_pos != std::string::npos) {
        rString.erase(end_pos, 1);
    }
}

/// @brief  Helper function to create a vector of strings from a string. ',' is used a delimiter.
/// @param rString
/// @return StringVectorType
inline StringVectorType CreateStringVector(const std::string& rString) {
    StringVectorType enum_strings;
    std::string cleanedString = GetStringWithoutWhiteSpaces(rString);
    RemoveQuESoList(cleanedString);
    std::istringstream stream(cleanedString);
    std::string token;
    while (std::getline(stream, token, ',')) {
        std::size_t pos = token.find('=');
        if (pos != std::string::npos) {
            QuESo_ERROR << "Keys must start at 0 and increment by 1. No custom values allowed.";
        } else {
            enum_strings.push_back(token);
        }
    }
    return enum_strings;
}

// Helper struct
template<typename... Types>
struct ValuesTypeList {
};

/// Forward declaration of KeyToWhatTypeTag template.
template<typename TDerived, typename TTypeList>
struct KeyToWhatTypeTag;

/// @brief This class acts as a type tag to static_assert, if a key is suiteable for certain objects.
/// @tparam TDerived - Curiously Recurring Template Pattern. This must be the derived class.
/// @tparam TTypeList - List of possible types for the underyling values. Example: You may create a
///         KeyToValues type tag that allows int and double as underlying values.
template<typename TDerived, typename... PossibleTypes>
struct KeyToWhatTypeTag<TDerived, ValuesTypeList<PossibleTypes...>> {

    // Delete copy/move constructor and assignment operator
    KeyToWhatTypeTag(const KeyToWhatTypeTag&) = delete;
    KeyToWhatTypeTag& operator=(const KeyToWhatTypeTag&) = delete;
    KeyToWhatTypeTag(KeyToWhatTypeTag&&) = delete;
    KeyToWhatTypeTag& operator=(KeyToWhatTypeTag&&) = delete;

    // Type traits to check the if a given type is a possible value type.
    template<typename TType>
    using is_valid_value_type = std::disjunction<std::is_same<TType, PossibleTypes>...>;

    template<typename TType>
    static constexpr bool is_valid_value_type_v = is_valid_value_type<TType>::value;

    /// @brief Returns the name of a given type.
    /// @tparam TType
    /// @return std::string
    template<typename TType>
    static std::string GetTypeName() {
        if constexpr(std::is_same_v<TType, TDerived>){ // To print name of derived class.
            return TDerived::msTypeName;
        }
        else if constexpr(KeyToWhatTypeTag::is_valid_value_type_v<TType>){ // To print name of values.
            return ValuesTypeList<PossibleTypes...>::msTypeNames[IndexOfType<TType, PossibleTypes...>::msIndex];
        }
        else {
            static_assert(KeyToWhatTypeTag::is_valid_value_type_v<TType>, "Invalid type provided to GetTypeName");
            return "";
        }
    }

private:
    // This is a pure static class.
    KeyToWhatTypeTag() {}
    ~KeyToWhatTypeTag() {}

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

} // End namespace detail
} // End namespace key
} // End namespace queso

namespace queso {
namespace key {

/// Base class of keys.
struct KeyBase {
    virtual ~KeyBase() = default;

    virtual IndexType Index() const = 0;
    virtual const std::string& Name() const = 0;
    virtual std::type_index VariableTypeIndex() const = 0;
    virtual std::type_index KeySetInfoTypeIndex() const = 0;
};

/// @brief Key class to define global keys.
/// @tparam TKeySetInfoType
/// @tparam TKeyToWhat
template<typename TKeySetInfoType, typename TKeyToWhat = typename TKeySetInfoType::KeySetToWhat>
struct Key : KeyBase {
    /// The following typedefs can be used for static type assertion.
    typedef TKeySetInfoType KeySetInfoType;
    typedef TKeyToWhat KeyToWhat;
    typedef typename KeySetInfoType::EnumType KeyValueType;

    /// @brief  Constructor.
    /// @param KeyValue
    Key(KeyValueType KeyValue) : mKeyValue(KeyValue) {
    }

    /// @brief Returns value/index of Key.
    /// @return IndexType
    IndexType Index() const override{
        return static_cast<IndexType>(mKeyValue);
    }

    /// @brief Returns name of Key.
    /// @return const std::string&
    const std::string& Name() const override {
        return KeySetInfoType::msEnumNames[static_cast<IndexType>(mKeyValue)];
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

private:
    /// Private member

    ///@brief Actual value of the key.
    const KeyValueType mKeyValue;
};

/// Base class to store and access information about a key set.
struct KeySetInfo {
    virtual ~KeySetInfo() = default;

    virtual const KeyBase* pGetKey(const std::string& rName) const = 0;
    virtual std::size_t GetNumberOfKeys() const noexcept = 0;
    virtual bool IsCorrectKeyType(const KeyBase& rKey) const noexcept = 0;
};

} // End namespace key
} // End namespace queso

#define QuESo_CREATE_VALUE_TYPE_LIST(ValuesTypeListName, ...) \
namespace key {\
    template<> \
    struct detail::ValuesTypeList<__VA_ARGS__> { \
        inline static const std::vector<std::string> msTypeNames = queso::key::detail::CreateStringVector(#__VA_ARGS__); \
    }; \
    typedef detail::ValuesTypeList<__VA_ARGS__> ValuesTypeListName; \
}\

#define QuESo_CREATE_KEY_SET_TO_OBJECT_TYPE_TAG(TypeTagName) \
namespace key {\
    struct TypeTagName : public detail::KeyToWhatTypeTag<TypeTagName, detail::ValuesTypeList<>> {\
        inline static const std::string msTypeName = #TypeTagName;\
    };\
}\

#define QuESo_CREATE_KEY_SET_TO_VALUE_TYPE_TAG(TypeTagName, ValueTypeList_) \
namespace key {\
    struct TypeTagName : public detail::KeyToWhatTypeTag<TypeTagName, ValueTypeList_> {\
        inline static const std::string msTypeName = #TypeTagName;\
    };\
}\

#define QuESo_KEY_LIST(...) __VA_ARGS__

#define QuESo_DEFINE_KEY_SET(KeySetName, KeySetToWhat_, ...) \
    typedef queso::key::detail::StringVectorType StringVectorType;\
    typedef queso::key::KeySetInfo KeySetInfoType;\
    typedef queso::key::KeyBase KeyBaseType;\
    namespace key {\
        /* KeySetName##KeyType##KeySetInfo allos to access key set information. */\
        struct KeySetName##KeySetToWhat_##KeySetInfo : public KeySetInfoType {\
            typedef queso::key::KeySetToWhat_ KeySetToWhat;\
            enum class EnumType {__VA_ARGS__};\
            /* Member function to access key information*/\
            const KeyBaseType* pGetKey(const std::string& rName) const override {\
                const auto it = msStringToKeyMap.find(rName);\
                if (it != msStringToKeyMap.end()) {\
                    return (it->second);\
                }\
                QuESo_ERROR << "Invalid Key name. Possible names are: " + StaticGetAllKeyNames();\
            }\
            std::size_t GetNumberOfKeys() const noexcept override {\
                return msEnumNames.size();\
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
            inline static const StringVectorType msEnumNames = queso::key::detail::CreateStringVector(#__VA_ARGS__);\
            static const std::unordered_map<std::string, const KeyBaseType* const> msStringToKeyMap; \
        };\
    }

#define QuESo_CREATE_KEY_SET_INFO(KeySetName, KeySetToWhat, KeyNames) \
    QuESo_DEFINE_KEY_SET(KeySetName, KeySetToWhat, KeyNames) \

#define QuESo_DEFINE_KEY_TO_OBJECT(KeySetName, KeyName, KeySetToWhat) \
namespace KeySetName {\
    inline const queso::key::Key<key::KeySetName##KeySetToWhat##KeySetInfo> KeyName = queso::key::Key<key::KeySetName##KeySetToWhat##KeySetInfo>(\
        key::KeySetName##KeySetToWhat##KeySetInfo::EnumType::KeyName);\
}\

#define QuESo_DEFINE_KEY_TO_VALUE(KeySetName, KeyName, KeySetToWhat, KeyToWhat) \
    namespace KeySetName {\
        static_assert(queso::key::KeySetToWhat::is_valid_value_type_v<KeyToWhat>, "Given KeyToWhat-type is invalid.");\
        inline const queso::key::Key<key::KeySetName##KeySetToWhat##KeySetInfo, KeyToWhat> KeyName  = queso::key::Key<key::KeySetName##KeySetToWhat##KeySetInfo, KeyToWhat>(\
            key::KeySetName##KeySetToWhat##KeySetInfo::EnumType::KeyName);\
    }\

namespace queso {
namespace key {
namespace detail {

typedef std::unordered_map<std::string, const queso::key::KeyBase* const> StringToKeyMapType;

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
            QuESo_ERROR << "Key '" + r_key_name + "' not found in key map. Check: REGISTER_KEYS_HPP.";
        }
    }
    return rKeyMap;
}

} // End namespace detail
} // End namespace key
} // End namespace queso

#define QuESo_KEY(Key) {Key.Name(), &Key} \

#define QuESo_REGISTER_KEY_SET(KeySet, KeySetToWhat, ...) \
    inline const queso::key::detail::StringToKeyMapType key::KeySet##KeySetToWhat##KeySetInfo::msStringToKeyMap =\
        queso::key::detail::InitializeStringToKeyMap<key::KeySet##KeySetToWhat##KeySetInfo>( {__VA_ARGS__} );\

#endif // End KEYS_INCLUDE_HPP