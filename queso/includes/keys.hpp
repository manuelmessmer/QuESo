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

//// External includes
#include "smart_enum/smartenum.hpp"

/// STL includes
#include <map> // Include map for ordered map


namespace queso {
namespace detail {

inline std::string RemoveWhiteSpace(const std::string& str) {
    std::string result;
    result.reserve(str.size());
    for (char ch : str) {
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

inline std::unordered_map<std::size_t, std::string> CreateEnumMap(const std::string& enumString) {
    std::unordered_map<std::size_t, std::string> enumMap;
    std::string cleanedString = RemoveWhiteSpace(enumString);
    RemoveQuESoList(cleanedString);
    std::istringstream stream(cleanedString);
    std::string token;
    std::size_t currentValue = 0;
    while (std::getline(stream, token, ',')) {
        std::size_t pos = token.find('=');
        if (pos != std::string::npos) {
            std::string name = token.substr(0, pos);
            std::size_t value = std::stoul(token.substr(pos + 1));
            enumMap[value] = name;
            currentValue = value + 1;
        } else {
            enumMap[currentValue] = token;
            currentValue++;
        }
    }
    return enumMap;
}

inline std::unordered_map<std::string, std::size_t> CreateReverseEnumMap(const std::unordered_map<std::size_t, std::string>& enumMap) {
    std::unordered_map<std::string, std::size_t> reverseEnumMap;
    for (const auto& pair : enumMap) {
        reverseEnumMap[pair.second] = pair.first;
    }
    return reverseEnumMap;
}

template <typename Type>
inline std::vector<Type> CreateEnumVector(const std::string& enumString) {
    std::vector<Type> enumVector;
    std::string cleanedString = RemoveWhiteSpace(enumString);
    RemoveQuESoList(cleanedString);
    std::istringstream stream(cleanedString);
    std::string token;
    while (std::getline(stream, token, ',')) {
        std::size_t pos = token.find('=');
        if (pos != std::string::npos) {
            std::string name = token.substr(0, pos);
            Type value = static_cast<Type>(std::stoul(token.substr(pos + 1)));
            enumVector.push_back(value);
        } else {
            enumVector.push_back(static_cast<Type>(enumVector.size()));
        }
    }
    return enumVector;
}

} // End namespace detail
} // End namespace queso

namespace queso {
namespace key {
    /// Type traits to distinguish between different key types.
    
    /// Key to list of DataSets. Each dataset (within one list) is accessed by Index not by key.
    struct List {
    };
    /// Key to subdictionary
    struct SubDict {
    };
    /// Key to dataset.
    struct DataSet {
    };

    /// Base class to access key information
    struct KeyInformation {
        virtual ~KeyInformation() = default;
        virtual const std::string& GetKeyName(std::size_t Index) const = 0;
        virtual std::size_t GetKeyValue(const std::string& rEnumName) const = 0;
        virtual std::size_t GetNumberOfKeys() const = 0;
        virtual const std::type_info& GetKeyTypeInfo() const = 0;
        virtual std::string GetAllKeyNames() const = 0;
        
    };
} // End namespace key
} // End namespace queso


#define QuESo_LIST(...) __VA_ARGS__

#define QuESo_CREATE_KEYS(KeyName, KeyType, Params, EnumType) \
    namespace key {\
        struct KeyName##KeyType##KeyInfo;\
        struct KeyName##KeyType {\
            using KeyToWhat = KeyType;\
            using KeyInfo = KeyName##KeyType##KeyInfo;\
        };\
        struct KeyName##KeyType##KeyInfo : public KeyInformation {\
            const std::string& GetKeyName(std::size_t Index) const override {\
                auto it = EnumNames.find(Index);\
                if (it != EnumNames.end()) {\
                    return it->second;\
                }\
                QuESo_ERROR << "Invalid index. Possible values are: " + StaticGetAllKeyNames();\
            }\
            std::size_t GetKeyValue(const std::string& rEnumName) const override {\
                auto it = EnumValues.find(rEnumName);\
                if (it != EnumValues.end()) {\
                    return it->second;\
                }\
                QuESo_ERROR << "Invalid enum name. Possible values are: " + StaticGetAllKeyNames();\
            }\
            std::size_t GetNumberOfKeys() const override {return Size;}\
            const std::type_info& GetKeyTypeInfo() const override {\
                return typeid(EnumType);\
            }\
            std::string GetAllKeyNames() const override {\
                return StaticGetAllKeyNames();\
            }\
            static std::string StaticGetAllKeyNames() {\
                std::string result;\
                std::map<std::size_t, std::string> ordered_names(EnumNames.begin(), EnumNames.end());\
                for (const auto& pair : ordered_names) {\
                    if (!result.empty()) {\
                        result += ", ";\
                    }\
                    result += pair.second;\
                }\
                return result.insert(0, 1, '[') + ']';\
            }\
            inline static std::unordered_map<std::size_t, std::string> EnumNames = queso::detail::CreateEnumMap(#Params);\
            inline static std::unordered_map<std::string, std::size_t> EnumValues = queso::detail::CreateReverseEnumMap(EnumNames);\
            inline static std::vector<EnumType> EnumList = queso::detail::CreateEnumVector<EnumType>(#Params);\
            inline static std::size_t Size = EnumList.size();\
        };\
        inline const std::string& KeyToString(EnumType value)\
        {\
            auto it = KeyName##KeyType##KeyInfo::EnumNames.find(static_cast<std::size_t>(value));\
            if (it != KeyName##KeyType##KeyInfo::EnumNames.end()) {\
                return it->second;\
            }\
            QuESo_ERROR << "Invalid enum value. Possible values are: " + KeyName##KeyType##KeyInfo::StaticGetAllKeyNames();\
        }\
        template<typename TType,\
                 typename = std::enable_if_t<std::is_same<TType, EnumType>::value>>\
        inline KeyName##KeyType GetKeyBaseType()\
        {\
            return KeyName##KeyType{};\
        }\
        inline std::ostream& operator<<(std::ostream& outStream, EnumType Value)\
        {\
            outStream << KeyToString(Value);\
            return outStream;\
        }\
        \
        template<typename TType,\
                 typename = std::enable_if_t<std::is_same<TType, KeyName##KeyType>::value>>\
        inline EnumType StringToKey(const std::string& rName)\
        {\
            const auto it = KeyName##KeyType##KeyInfo::EnumValues.find(rName);\
            if (it != KeyName##KeyType##KeyInfo::EnumValues.end()) {\
                return static_cast<EnumType>(it->second);\
            }\
            QuESo_ERROR << "Invalid enum name. Possible values are: " + KeyName##KeyType##KeyInfo::StaticGetAllKeyNames();\
        }\
        inline bool IsCorrectType(Unique<KeyInformation>& pKeyInformation, EnumType Value) {\
            return pKeyInformation->GetKeyTypeInfo() == typeid(Value);\
        }\
    }

#define QuESo_REGISTER_KEYS_1(KeyName, KeyType, Params) \
    struct KeyName {\
        enum KeyTo##KeyType {Params};\
    };\
    QuESo_CREATE_KEYS(KeyName, KeyType, QuESo_LIST(Params), KeyName::KeyTo##KeyType ) \

#define QuESo_REGISTER_KEYS_2(KeyName, KeyType1, Params1, KeyType2, Params2) \
    namespace key {\
    namespace detail {\
    enum class CheckDuplicatedValuesOf##KeyName {Params1, Params2};\
    }\
    }\
    struct KeyName {\
        enum KeyTo##KeyType1 {Params1};\
        enum KeyTo##KeyType2 {Params2};\
    };\
    QuESo_CREATE_KEYS(KeyName, KeyType1, QuESo_LIST(Params1), KeyName::KeyTo##KeyType1)\
    QuESo_CREATE_KEYS(KeyName, KeyType2, QuESo_LIST(Params2), KeyName::KeyTo##KeyType2)\

#define QuESo_REGISTER_KEYS_3(KeyName, KeyType1, Params1, KeyType2, Params2, KeyType3, Params3)\
    namespace key {\
    namespace detail {\
    enum class CheckDuplicatedValuesOf##KeyName {Params1, Params2, Params3};\
    }\
    }\
    struct KeyName {\
        enum KeyTo##KeyType1 {Params1};\
        enum KeyTo##KeyType2 {Params2};\
        enum KeyTo##KeyType2 {Params3};\
    };\
    QuESo_CREATE_KEYS(KeyName, KeyType1, QuESo_LIST(Params1), KeyName::KeyTo##KeyType1)\
    QuESo_CREATE_KEYS(KeyName, KeyType2, QuESo_LIST(Params2), KeyName::KeyTo##KeyType2)\
    QuESo_CREATE_KEYS(KeyName, KeyType3, QuESo_LIST(Params3), KeyName::KeyTo##KeyType3)\

#endif // End KEYS_INCLUDE_HPP