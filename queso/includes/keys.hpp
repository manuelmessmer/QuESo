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
#include <optional>


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
        virtual const std::string& EnumName(std::size_t Index) = 0;
        virtual std::size_t EnumValue(const std::string& rEnumName) = 0;
        virtual std::size_t EnumSize() = 0;
    };
} // End namesapce key
} // End namespace queso


#define QuESo_LIST(...) __VA_ARGS__

#define QuESo_CREATE_KEYS(KeyName, KeyType, Params) \
    namespace key {\
    struct KeyName##KeyType : public KeyInformation {\
        using type = KeyType;\
        enum BaseEnum {Params};\
        const std::string& EnumName(std::size_t Index) override {\
            auto it = EnumNames.find(Index);\
            if (it != EnumNames.end()) {\
                return it->second;\
            }\
            QuESo_ERROR << "Invalid index. Possible values are: " + GetAllEnumNames();\
        }\
        std::size_t EnumValue(const std::string& rEnumName) override {\
            auto it = EnumValues.find(rEnumName);\
            if (it != EnumValues.end()) {\
                return it->second;\
            }\
            QuESo_ERROR << "Invalid enum name. Possible values are: " + GetAllEnumNames();\
        }\
        std::size_t EnumSize() override {return Size;}\
        inline static std::unordered_map<std::size_t, std::string> EnumNames = queso::detail::CreateEnumMap(#Params);\
        inline static std::unordered_map<std::string, std::size_t> EnumValues = queso::detail::CreateReverseEnumMap(EnumNames);\
        inline static std::vector<BaseEnum> EnumList = queso::detail::CreateEnumVector<BaseEnum>(#Params);\
        inline static std::size_t Size = EnumList.size();\
        static std::string GetAllEnumNames() {\
            std::string result;\
            for (const auto& pair : EnumNames) {\
                if (!result.empty()) {\
                    result += ", ";\
                }\
                result += pair.second;\
            }\
            return result;\
        }\
    };\
    inline const std::string& KeyToString(KeyName##KeyType::BaseEnum value) \
    { \
        auto it = KeyName##KeyType::EnumNames.find(static_cast<std::size_t>(value));\
        if (it != KeyName##KeyType::EnumNames.end()) {\
            return it->second;\
        }\
        QuESo_ERROR << "Invalid enum value. Possible values are: " + KeyName##KeyType::GetAllEnumNames();\
    }\
    inline KeyName##KeyType::type GetType(KeyName##KeyType::BaseEnum A)\
    {\
        return KeyName##KeyType::type{};\
    }\
    inline std::ostream& operator<<(std::ostream& outStream, KeyName##KeyType::BaseEnum value)\
    {\
        outStream << KeyToString(value);\
        return outStream;\
    }\
    \
    template<typename TType,\
             typename = std::enable_if_t<std::is_same<TType, KeyName##KeyType::BaseEnum>::value>>\
    inline const KeyName##KeyType::BaseEnum StringToKey(const std::string& name)\
    {\
        auto it = KeyName##KeyType::EnumValues.find(name);\
        if (it != KeyName##KeyType::EnumValues.end()) {\
            return static_cast<KeyName##KeyType::BaseEnum>(it->second);\
        }\
        QuESo_ERROR << "Invalid enum name. Possible values are: " + KeyName##KeyType::GetAllEnumNames();\
    }\
    } // End namespace key

#define QuESo_REGISTER_KEYS_1(KeyName, KeyType, Params) \
    QuESo_CREATE_KEYS(KeyName, KeyType, QuESo_LIST(Params)) \
    typedef key::KeyName##KeyType KeyName;

#define QuESo_REGISTER_KEYS_2(KeyName, KeyType1, Params1, KeyType2, Params2) \
    namespace detail {\
    enum class CheckDuplicatedValuesOf##KeyName {Params1, Params2}; \
    }\
    QuESo_CREATE_KEYS(KeyName, KeyType1, QuESo_LIST(Params1)) \
    QuESo_CREATE_KEYS(KeyName, KeyType2, QuESo_LIST(Params2)) \
    struct KeyName : public key::KeyName##KeyType1, key::KeyName##KeyType2 {};

#define QuESo_REGISTER_KEYS_3(KeyName, KeyType1, Params1, KeyType2, Params2, KeyType3, Params3) \
    namespace detail {\
    enum class CheckDuplicatedValuesOf##KeyName {Params1, Params2, Params3}; \
    }\
    QuESo_CREATE_KEYS(KeyName, KeyType1, QuESo_LIST(Params1)) \
    QuESo_CREATE_KEYS(KeyName, KeyType2, QuESo_LIST(Params2)) \
    QuESo_CREATE_KEYS(KeyName, KeyType3, QuESo_LIST(Params3)) \
    struct KeyName : public key::KeyName##KeyType1, key::KeyName##KeyType2, key::KeyName##KeyType3 {};

#endif // KEYS_INCLUDE