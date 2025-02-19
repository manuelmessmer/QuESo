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
#include <map> 
#include <typeindex>

namespace queso {
namespace key {
namespace detail {

typedef std::vector<std::string> StringVectorType;

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

inline void RemoveQuESoList(std::string& rString) {
    std::size_t start_pos = rString.find("QuESo_LIST(");
    if (start_pos != std::string::npos) {
        rString.erase(start_pos, std::string("QuESo_LIST(").length());
    }
    std::size_t end_pos = rString.find(")", start_pos);
    if (end_pos != std::string::npos) {
        rString.erase(end_pos, 1);
    }
}

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

template<typename TEnumType>
inline std::unordered_map<std::string, TEnumType> CreateStringToEnumMap(const StringVectorType& rStringVector) {
    std::unordered_map<std::string, TEnumType> reverse_enum_map;
    for (std::size_t i = 0; i < rStringVector.size(); ++i) {
        reverse_enum_map[rStringVector[i]] = static_cast<TEnumType>(i);
    }
    return reverse_enum_map;
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
    /// Key to dataset.
    struct KeyToValue {
    };

    struct KeyBase {
        virtual IndexType Index() const = 0;
        virtual const std::string& Name() const = 0;
        virtual std::type_index VariableTypeInfo() const = 0;
    };

    template<typename TKeyInfoType, typename TKeyToWhat = typename TKeyInfoType::KeyToWhat>
    struct Key : KeyBase {
        typedef TKeyInfoType KeyInfoType;
        typedef TKeyToWhat KeyToWhat;
        typedef typename KeyInfoType::EnumType KeyValueType;

        Key(KeyValueType KeyValue) : mKeyValue(KeyValue) {
        }

        Key(const Key &rOther) = delete;

        IndexType Index() const override{
            return static_cast<IndexType>(mKeyValue);
        }
        
        const std::string& Name() const override {
            return KeyInfoType::msEnumNames[static_cast<IndexType>(mKeyValue)];
        }

        std::type_index VariableTypeInfo() const override {
            return std::type_index(typeid(KeyToWhat));
        }

        KeyValueType Value(){
            return mKeyValue;
        }
    private:
        const KeyValueType mKeyValue;
    };

    /// Base class to store and access key information.
    struct KeyInformation {
        virtual ~KeyInformation() = default;
        virtual const KeyBase* GetKey(const std::string& rName) const = 0;
        virtual std::size_t GetNumberOfKeys() const noexcept = 0;        
    };
} // End namespace key
} // End namespace queso

#define QuESo_LIST(...) __VA_ARGS__

#define QuESo_CREATE_KEY_INFO(KeySetName, KeyToWhat_, KeyNames) \
    typedef queso::key::detail::StringVectorType StringVectorType;\
    typedef queso::key::KeyInformation KeyInformationType;\
    typedef queso::key::KeyBase KeyBaseType;\
    namespace key {\
        /* KeySetName##KeyType##KeyInfo allos to access enum/key information. */\
        struct KeySetName##KeyToWhat_##KeyInfo : public KeyInformationType {\
            typedef queso::key::KeyToWhat_ KeyToWhat;\
            enum class EnumType {KeyNames};\
            /* Member function to access enum/key information*/\
            const KeyBaseType* GetKey(const std::string& rName) const override {\
                const auto it = msStringToKeyMap.find(rName);\
                if (it != msStringToKeyMap.end()) {\
                    return (it->second);\
                }\
                QuESo_ERROR << "Invalid Key name. Possible names are: " + StaticGetAllKeyNames();\
            }\
            std::size_t GetNumberOfKeys() const noexcept override {\
                return msEnumNames.size();\
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
            /* Static members that contain enum/key information, e.g., to map enum to string, etc. */\
            inline static const StringVectorType msEnumNames = queso::key::detail::CreateStringVector(#KeyNames);\
            static const std::unordered_map<std::string, KeyBaseType*> msStringToKeyMap; \
        };\
    }

#define QuESo_DEFINE_KEY_SET(KeySetName, KeySetToWhat, KeyNames) \
    QuESo_CREATE_KEY_INFO(KeySetName, KeySetToWhat, QuESo_LIST(KeyNames)) \

#define QuESo_DEFINE_KEY_TO_VALUE(KeySetName, KeyName, KeyToWhat) \
    namespace KeySetName {\
        inline const queso::key::Key<key::KeySetName##KeyToValue##KeyInfo, KeyToWhat> KeyName  = queso::key::Key<key::KeySetName##KeyToValue##KeyInfo, KeyToWhat>(\
            key::KeySetName##KeyToValue##KeyInfo::EnumType::KeyName);\
    }\

#define QuESo_DEFINE_KEY_TO_SUBDICT(KeySetName, KeyName) \
    namespace KeySetName {\
        extern queso::key::Key<key::KeySetName##KeyToSubDict##KeyInfo> KeyName;\
    }\

#define QuESo_DEFINE_KEY_TO_LIST(KeySetName, KeyName) \
    namespace KeySetName {\
        extern queso::key::Key<key::KeySetName##KeyToList##KeyInfo> KeyName;\
    }\

#define QuESo_CREATE_KEY(KeySetName, KeyName) \
    namespace KeySetName {\
        queso::key::Key<decltype(KeyName)::KeyInfoType, decltype(KeyName)::KeyToWhat> KeyName(decltype(KeyName)::KeyInfoType::EnumType::KeyName);\
    }\

#define QuESo_KEY(Key) {Key.Name(), &Key} \

#define QuESo_REGISTER_KEY_SET(KeySet, KeySetToWhat, ...) \
    const std::unordered_map<std::string, queso::key::KeyBase*> key::KeySet##KeySetToWhat##KeyInfo::msStringToKeyMap = { __VA_ARGS__ };\



#endif // End KEYS_INCLUDE_HPP