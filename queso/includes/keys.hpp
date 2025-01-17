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
#include <unordered_map>
#include <map> 

namespace queso {
namespace key {
namespace detail {

typedef std::vector<std::string> StringVectorType;
typedef std::unordered_map<std::string, std::size_t> StringToIndexMapType;

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

inline StringToIndexMapType CreateStringToIndexMap(const StringVectorType& rStringVector) {
    StringToIndexMapType reverse_enum_map;
    for (std::size_t i = 0; i < rStringVector.size(); ++i) {
        reverse_enum_map[rStringVector[i]] = i;
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
    struct List {
    };
    /// Key to subdictionary
    struct SubDict {
    };
    /// Key to dataset.
    struct DataSet {
    };

    /// Base class to store and access key information.
    struct KeyInformation {
        virtual ~KeyInformation() = default;
        virtual const std::string& GetKeyName(std::size_t Index) const = 0;
        virtual std::size_t GetKeyValue(const std::string& rEnumName) const = 0;
        virtual std::size_t GetNumberOfKeys() const noexcept = 0;
        virtual const std::type_info& GetKeyTypeInfo() const = 0;
        virtual std::string GetAllKeyNames() const = 0;
        
    };
} // End namespace key
} // End namespace queso


#define QuESo_LIST(...) __VA_ARGS__

#define QuESo_CREATE_KEY_INFO(KeySetName, KeyType, KeyNames, EnumType) \
    typedef queso::key::detail::StringVectorType StringVectorType;\
    typedef queso::key::detail::StringToIndexMapType StringToIndexMapType;\
    namespace key {\
        struct KeySetName##KeyType##KeyInfo; /* Forward declaration */\
        struct KeySetName##KeyType {\
            using KeyToWhat = KeyType;\
            using KeyInfo = KeySetName##KeyType##KeyInfo;\
        };\
        /* KeySetName##KeyType##KeyInfo allos to access enum/key information. */\
        struct KeySetName##KeyType##KeyInfo : public KeyInformation {\
            /* Member function to access enum/key information*/\
            const std::string& GetKeyName(std::size_t Index) const override {\
                if( Index >= 0 && Index < msSize ) {\
                    return msEnumNames[Index];\
                }\
                QuESo_ERROR << "Invalid index. Possible values are: " + StaticGetAllKeyNames();\
            }\
            std::size_t GetKeyValue(const std::string& rEnumName) const override {\
                const auto it = msEnumValues.find(rEnumName);\
                if (it != msEnumValues.end()) {\
                    return it->second;\
                }\
                QuESo_ERROR << "Invalid Key name. Possible names are: " + StaticGetAllKeyNames();\
            }\
            std::size_t GetNumberOfKeys() const noexcept override {\
                return msSize;\
            }\
            const std::type_info& GetKeyTypeInfo() const override {\
                return typeid(EnumType);\
            }\
            std::string GetAllKeyNames() const override {\
                return StaticGetAllKeyNames();\
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
            inline static StringVectorType msEnumNames = queso::key::detail::CreateStringVector(#KeyNames);\
            inline static StringToIndexMapType msEnumValues = queso::key::detail::CreateStringToIndexMap(msEnumNames);\
            inline static std::size_t msSize = msEnumNames.size();\
        };\
        inline const std::string& KeyToString(EnumType value) {\
            if( static_cast<std::size_t>(value) >= 0\
                    && static_cast<std::size_t>(value) < KeySetName##KeyType##KeyInfo::msEnumNames.size() ) {\
                return KeySetName##KeyType##KeyInfo::msEnumNames[static_cast<std::size_t>(value)];\
            }\
            QuESo_ERROR << "Invalid enum value. Possible values are: " + KeySetName##KeyType##KeyInfo::StaticGetAllKeyNames();\
        }\
        template<typename TType,\
                 typename = std::enable_if_t<std::is_same<TType, KeySetName##KeyType>::value>>\
        inline EnumType StringToKey(const std::string& rName) {\
            const auto it = KeySetName##KeyType##KeyInfo::msEnumValues.find(rName);\
            if (it != KeySetName##KeyType##KeyInfo::msEnumValues.end()) {\
                return static_cast<EnumType>(it->second);\
            }\
            QuESo_ERROR << "Invalid enum name. Possible values are: " + KeySetName##KeyType##KeyInfo::StaticGetAllKeyNames();\
        }\
        inline bool IsCorrectType(Unique<KeyInformation>& pKeyInformation, EnumType Value) {\
            return pKeyInformation->GetKeyTypeInfo() == typeid(Value);\
        }\
        template<typename TType,\
                 typename = std::enable_if_t<std::is_same<TType, EnumType>::value>>\
        inline KeySetName##KeyType GetKeyBaseType() noexcept {\
            return KeySetName##KeyType{};\
        }\
        inline std::ostream& operator<<(std::ostream& outStream, EnumType Value) {\
            outStream << KeyToString(Value);\
            return outStream;\
        }\
    }

/**
 * @brief Macro to register a key set for a given key set name and key type.
 *
 * This macro defines an enum within a struct for the given key set name and key type
 * and creates the respective KeyInformation (QuESo_CREATE_KEY_INFO macro), which allows
 * to get more information about the enum/key set, e.g., KeyInformation allows to map from key to string.
 *
 * @param KeySetName The name of the key set.
 * @param KeyType The type of the key (possible options List, SubDict, DataSet).
 * @param KeyNames The names of the actual keys, specified as a QuESo_LIST (QuESo_LIST macro).
 *
 * Example usage:
 * @code
 * QuESo_REGISTER_KEYS_1(KeySetName, List, QuESo_LIST(KeyName1 = 0, KeyName2, KeyName3 = 5) )
 * @endcode
 */
#define QuESo_REGISTER_KEY_SET_1(KeySetName, KeyType, KeyNames) \
    struct KeySetName {\
        enum KeyTo##KeyType {KeyNames};\
    };\
    QuESo_CREATE_KEY_INFO(KeySetName, KeyType, QuESo_LIST(KeyNames), KeySetName::KeyTo##KeyType ) \

/**
 * @brief Macro to register a key set for a given key set name and two key types.
 * 
 * This macro defines two enums within a struct for the given key set name and key types
 * and creates the respective KeyInformation (QuESo_CREATE_KEY_INFO macro), which allows
 * to get more information about the enum/key set, e.g., KeyInformation allows to map 
 * from key to string. 
 * For each key type a different enum is created. This allows to distinguish between keys 
 * that access e.g. List or DataSets. 
 *
 * @param KeySetName The name of the key set.
 * @param KeyType1 The type of the first key collection (possible options List, SubDict, DataSet).
 * @param KeyNames1 The names of the actual keys for the first key collection, specified as a QuESo_LIST (QuESo_LIST macro).
 * @param KeyType2 The type of the second key (possible options List, SubDict, DataSet).
 * @param KeyNames2 The names of the actual keys for the second key type, specified as a QuESo_LIST (QuESo_LIST macro).
 *
 * Example usage:
 * @code
 * QuESo_REGISTER_KEYS_2(KeySetName, List, QuESo_LIST(KeyName1 = 0, KeyName2), SubDict, QuESo_LIST(KeyName3, KeyName4 = 5) )
 * @endcode
 */
#define QuESo_REGISTER_KEY_SET_2(KeySetName, KeyType1, KeyNames1, KeyType2, KeyNames2) \
    namespace key {\
    namespace detail {\
        enum class CheckDuplicatedValuesOf##KeySetName {KeyNames1, KeyNames2};\
    }\
    }\
    struct KeySetName {\
        enum KeyTo##KeyType1 {KeyNames1};\
        enum KeyTo##KeyType2 {KeyNames2};\
    };\
    QuESo_CREATE_KEY_INFO(KeySetName, KeyType1, QuESo_LIST(KeyNames1), KeySetName::KeyTo##KeyType1)\
    QuESo_CREATE_KEY_INFO(KeySetName, KeyType2, QuESo_LIST(KeyNames2), KeySetName::KeyTo##KeyType2)\

/**
 * @brief Macro to register a key set for a given key set name and three key types.
 * 
 * This macro defines three enums within a struct for the given key set name and key types
 * and creates the respective KeyInformation (QuESo_CREATE_KEY_INFO macro), which allows
 * to get more information about the enum/key set, e.g., KeyInformation allows to map 
 * from key to string. 
 * For each key type a different enum is created. This allows to distinguish between keys 
 * that access e.g. a List or a DataSet. 
 *
 * @param KeySetName The name of the key set.
 * @param KeyType1 The type of the first key collection (possible options List, SubDict, DataSet).
 * @param KeyNames1 The names of the actual keys for the first key collection, specified as a QuESo_LIST (QuESo_LIST macro).
 * @param KeyType2 The type of the second key (possible options List, SubDict, DataSet).
 * @param KeyNames2 The names of the actual keys for the second key type, specified as a QuESo_LIST (QuESo_LIST macro).
 * @param KeyType3 The type of the third key (possible options List, SubDict, DataSet).
 * @param KeyNames3 The names of the actual keys for the third key type, specified as a QuESo_LIST (QuESo_LIST macro).
 *
 * Example usage:
 * @code
 * QuESo_REGISTER_KEYS_3(KeySetName, List, QuESo_LIST(KeyName1 = 0, KeyName2), SubDict, QuESo_LIST(KeyName3, KeyName4 = 5), DataSet, QuESo_LIST(KeyName5, KeyName6 = 10) )
 * @endcode
 */
#define QuESo_REGISTER_KEY_SET_3(KeySetName, KeyType1, KeyNames1, KeyType2, KeyNames2, KeyType3, KeyNames3)\
    namespace key {\
    namespace detail {\
        enum class CheckDuplicatedValuesOf##KeySetName {KeyNames1, KeyNames2, KeyNames3};\
    }\
    }\
    struct KeySetName {\
        enum KeyTo##KeyType1 {KeyNames1};\
        enum KeyTo##KeyType2 {KeyNames2};\
        enum KeyTo##KeyType3 {KeyNames3};\
    };\
    QuESo_CREATE_KEY_INFO(KeySetName, KeyType1, QuESo_LIST(KeyNames1), KeySetName::KeyTo##KeyType1)\
    QuESo_CREATE_KEY_INFO(KeySetName, KeyType2, QuESo_LIST(KeyNames2), KeySetName::KeyTo##KeyType2)\
    QuESo_CREATE_KEY_INFO(KeySetName, KeyType3, QuESo_LIST(KeyNames3), KeySetName::KeyTo##KeyType3)\

#endif // End KEYS_INCLUDE_HPP