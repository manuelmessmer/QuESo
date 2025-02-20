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

} // End namespace detail
} // End namespace key
} // End namespace queso

namespace queso {
namespace key {
    /// Type traits to distinguish between different key types(e.g., List, SubDict, DataSet).

    /// Key to list of DataSets. Each dataset (within one list) is accessed by Index not by key.
    struct KeyToList {
    };
    /// Key to subdictionary
    struct KeyToSubDict {
    };
    /// Key to value.
    struct KeyToValue {
    };

    /// Base class of keys.
    struct KeyBase {
        virtual IndexType Index() const = 0;
        virtual const std::string& Name() const = 0;
        virtual std::type_index VariableTypeIndex() const = 0;
        virtual std::type_index KeySetInfoTypeIndex() const = 0;
    };

    /// @brief Key class to define global keys.
    /// @tparam TKeySetInfoType
    /// @tparam TKeyToWhat
    template<typename TKeySetInfoType, typename TKeyToWhat = typename TKeySetInfoType::KeyToWhat>
    struct Key : KeyBase {
        /// The following typedefs can be used for static type checks.
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

#define QuESo_KEY_LIST(...) __VA_ARGS__

#define QuESo_CREATE_KEY_SET_INFO(KeySetName, KeyToWhat_, KeyNames) \
    typedef queso::key::detail::StringVectorType StringVectorType;\
    typedef queso::key::KeySetInfo KeySetInfoType;\
    typedef queso::key::KeyBase KeyBaseType;\
    namespace key {\
        /* KeySetName##KeyType##KeySetInfo allos to access key set information. */\
        struct KeySetName##KeyToWhat_##KeySetInfo : public KeySetInfoType {\
            typedef queso::key::KeyToWhat_ KeyToWhat;\
            enum class EnumType {KeyNames};\
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
                return (std::type_index(typeid(KeySetName##KeyToWhat_##KeySetInfo)) == rKey.KeySetInfoTypeIndex() );\
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
            inline static const StringVectorType msEnumNames = queso::key::detail::CreateStringVector(#KeyNames);\
            static const std::unordered_map<std::string, const KeyBaseType* const> msStringToKeyMap; \
        };\
    }

#define QuESo_DEFINE_KEY_SET(KeySetName, KeySetToWhat, KeyNames) \
    QuESo_CREATE_KEY_SET_INFO(KeySetName, KeySetToWhat, QuESo_KEY_LIST(KeyNames)) \

#define QuESo_DEFINE_KEY_TO_VALUE(KeySetName, KeyName, KeyToWhat) \
    namespace KeySetName {\
        inline const queso::key::Key<key::KeySetName##KeyToValue##KeySetInfo, KeyToWhat> KeyName  = queso::key::Key<key::KeySetName##KeyToValue##KeySetInfo, KeyToWhat>(\
            key::KeySetName##KeyToValue##KeySetInfo::EnumType::KeyName);\
    }\

#define QuESo_DEFINE_KEY_TO_SUBDICT(KeySetName, KeyName) \
    namespace KeySetName {\
        inline const queso::key::Key<key::KeySetName##KeyToSubDict##KeySetInfo> KeyName = queso::key::Key<key::KeySetName##KeyToSubDict##KeySetInfo>(\
            key::KeySetName##KeyToSubDict##KeySetInfo::EnumType::KeyName);\
    }\

#define QuESo_DEFINE_KEY_TO_LIST(KeySetName, KeyName) \
    namespace KeySetName {\
        inline const queso::key::Key<key::KeySetName##KeyToList##KeySetInfo> KeyName = queso::key::Key<key::KeySetName##KeyToList##KeySetInfo>(\
            key::KeySetName##KeyToList##KeySetInfo::EnumType::KeyName);\
    }\

#define QuESo_KEY(Key) {Key.Name(), &Key} \

namespace queso {
namespace key {
namespace detail {
    typedef std::unordered_map<std::string, const queso::key::KeyBase* const> StringToKeyMapType;

    /// @brief Helper functions to check if all keys are correctly defined.
    /// @tparam TKeySetInfoType
    /// @param rKeyMap
    /// @return StringToKeyMapType
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
            if( it->first != r_key_name) {
                QuESo_ERROR << "Key '" + r_key_name + "' not found in key map. Check: REGISTER_KEYS_HPP.";
            }
        }
        return rKeyMap;
    }

} // End namespace detail
} // End namespace key
} // End namespace queso

#define QuESo_REGISTER_KEY_SET(KeySet, KeySetToWhat, ...) \
    inline const queso::key::detail::StringToKeyMapType key::KeySet##KeySetToWhat##KeySetInfo::msStringToKeyMap =\
        queso::key::detail::InitializeStringToKeyMap<key::KeySet##KeySetToWhat##KeySetInfo>( {__VA_ARGS__} );\

#endif // End KEYS_INCLUDE_HPP