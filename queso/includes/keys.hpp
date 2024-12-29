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

inline bool ReplaceString(std::string& rString, const std::string& rFrom, const std::string& rTo) {
    std::size_t start_pos = rString.find(rFrom);
    if(start_pos == std::string::npos)
        return false;
    rString.replace(start_pos, rFrom.length(), rTo);
    return true;
}

inline std::string PrepareKeyStrings(std::string enumValuesString)
{
    if( ReplaceString(enumValuesString, "QuESo_LIST(", "") 
        || ReplaceString(enumValuesString, ")", "") ) {
            QuESo_ERROR << "Strings could not be prepared succesfully.\n";
    }

    return enumValuesString;
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
}
} // End namespace queso


#define QuESo_LIST(...) __VA_ARGS__

#define QuESo_CREATE_KEYS(KeyName, KeyType, Params) \
    namespace key {\
    struct KeyName##KeyType : public KeyInformation {\
        using type = KeyType;\
        enum BaseEnum {Params};\
        const std::string& EnumName(std::size_t Index) override {return EnumNames.at(Index);}\
        std::size_t EnumValue(const std::string& rEnumName) override {return EnumValues.at(rEnumName);}\
        std::size_t EnumSize() override {return Size;}\
        inline static std::unordered_map<std::size_t, std::string> EnumNames = smart_enum::makeEnumNameMap(queso::detail::PrepareKeyStrings(#Params));\
        inline static std::vector<BaseEnum> EnumList = smart_enum::makeEnumList<BaseEnum>(queso::detail::PrepareKeyStrings(#Params));\
        inline static std::unordered_map<std::string, std::size_t> EnumValues = smart_enum::makeEnumValuesMap(queso::detail::PrepareKeyStrings(#Params));\
        inline static std::size_t Size = EnumList.size();\
    };\
    inline const std::string& KeyToString(KeyName##KeyType::BaseEnum value) \
    { \
        return KeyName##KeyType::EnumNames.at(static_cast<std::size_t>(value));\
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
        return static_cast<KeyName##KeyType::BaseEnum>(KeyName##KeyType::EnumValues.at(name));\
    }\
    }

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

#endif // KEYS_INCLUDE_HPP