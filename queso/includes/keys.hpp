// //   ____        ______  _____
// //  / __ \      |  ____|/ ____|
// // | |  | |_   _| |__  | (___   ___
// // | |  | | | | |  __|  \___ \ / _ \'
// // | |__| | |_| | |____ ____) | (_) |
// //  \___\_\\__,_|______|_____/ \___/
// //         Quadrature for Embedded Solids
// //
// //  License:    BSD 4-Clause License
// //              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
// //
// //  Authors:    Manuel Messmer
// //

// #ifndef KEYS_INCLUDE_HPP
// #define KEYS_INCLUDE_HPP

// /// STL includes
// #include <unordered_map>
// #include <map> 

// namespace queso {
// namespace key {
// namespace detail {

// typedef std::unordered_map<std::size_t, std::string> EnumMapType;
// typedef std::unordered_map<std::string, std::size_t> ReverseEnumMapType;

// inline std::string GetStringWithoutWithSpace(const std::string& rString) {
//     std::string result;
//     result.reserve(rString.size());
//     for (char ch : rString) {
//         if (!std::isspace(static_cast<unsigned char>(ch))) {
//             result.push_back(ch);
//         }
//     }
//     return result;
// }

// inline void RemoveQuESoList(std::string& rString) {
//     std::size_t start_pos = rString.find("QuESo_LIST(");
//     if (start_pos != std::string::npos) {
//         rString.erase(start_pos, std::string("QuESo_LIST(").length());
//     }
//     std::size_t end_pos = rString.find(")", start_pos);
//     if (end_pos != std::string::npos) {
//         rString.erase(end_pos, 1);
//     }
// }

// inline EnumMapType CreateEnumMap(const std::string& rString) {
//     EnumMapType enum_map;
//     std::string cleanedString = GetStringWithoutWithSpace(rString);
//     RemoveQuESoList(cleanedString);
//     std::istringstream stream(cleanedString);
//     std::string token;
//     std::size_t currentValue = 0;
//     while (std::getline(stream, token, ',')) {
//         std::size_t pos = token.find('=');
//         if (pos != std::string::npos) {
//             std::string name = token.substr(0, pos);
//             std::size_t value = std::stoul(token.substr(pos + 1));
//             enum_map[value] = name;
//             currentValue = value + 1;
//         } else {
//             enum_map[currentValue] = token;
//             currentValue++;
//         }
//     }
//     return enum_map;
// }

// inline ReverseEnumMapType CreateReverseEnumMap(const EnumMapType& rEnumMap) {
//     ReverseEnumMapType reverse_enum_map;
//     for (const auto& r_pair : rEnumMap) {
//         reverse_enum_map[r_pair.second] = r_pair.first;
//     }
//     return reverse_enum_map;
// }

// } // End namespace detail
// } // End namespace key
// } // End namespace queso

// namespace queso {
// namespace key {
//     /// Type traits to distinguish between different key types(e.g., List, SubDict, DataSet).
    
//     /// Key to list of DataSets. Each dataset (within one list) is accessed by Index not by key.
//     struct List {
//     };
//     /// Key to subdictionary
//     struct SubDict {
//     };
//     /// Key to dataset.
//     struct DataSet {
//     };

//     /// Base class to store and access key information.
//     struct KeyInformation {
//         virtual ~KeyInformation() = default;
//         virtual const std::string& GetKeyName(std::size_t Index) const = 0;
//         virtual std::size_t GetKeyValue(const std::string& rEnumName) const = 0;
//         virtual std::size_t GetNumberOfKeys() const noexcept = 0;
//         virtual const std::type_info& GetKeyTypeInfo() const = 0;
//         virtual std::string GetAllKeyNames() const = 0;
        
//     };
// } // End namespace key
// } // End namespace queso


// #define QuESo_LIST(...) __VA_ARGS__

// #define QuESo_CREATE_KEY_INFO(KesSetName, KeyType, KeyNames, EnumType) \
//     typedef queso::key::detail::EnumMapType EnumMapType;\
//     typedef queso::key::detail::ReverseEnumMapType ReverseEnumMapType;\
//     namespace key {\
//         struct KesSetName##KeyType##KeyInfo; /* Forward declaration */\
//         struct KesSetName##KeyType {\
//             using KeyToWhat = KeyType;\
//             using KeyInfo = KesSetName##KeyType##KeyInfo;\
//         };\
//         /* KesSetName##KeyType##KeyInfo allos to access enum/key information. */\
//         struct KesSetName##KeyType##KeyInfo : public KeyInformation {\
//             /* Member function to access enum/key information*/\
//             const std::string& GetKeyName(std::size_t Index) const override {\
//                 const auto it = msEnumNames.find(Index);\
//                 if (it != msEnumNames.end()) {\
//                     return it->second;\
//                 }\
//                 QuESo_ERROR << "Invalid index. Possible values are: " + StaticGetAllKeyNames();\
//             }\
//             std::size_t GetKeyValue(const std::string& rEnumName) const override {\
//                 const auto it = msEnumValues.find(rEnumName);\
//                 if (it != msEnumValues.end()) {\
//                     return it->second;\
//                 }\
//                 QuESo_ERROR << "Invalid Key name. Possible names are: " + StaticGetAllKeyNames();\
//             }\
//             std::size_t GetNumberOfKeys() const noexcept override {\
//                 return msSize;\
//             }\
//             const std::type_info& GetKeyTypeInfo() const override {\
//                 return typeid(EnumType);\
//             }\
//             std::string GetAllKeyNames() const override {\
//                 return StaticGetAllKeyNames();\
//             }\
//             static std::string StaticGetAllKeyNames() {\
//                 std::string result;\
//                 std::map<std::size_t, std::string> ordered_names(msEnumNames.begin(), msEnumNames.end());\
//                 for (const auto& pair : ordered_names) {\
//                     if (!result.empty()) {\
//                         result += "', '";\
//                     }\
//                     result += pair.second;\
//                 }\
//                 return result.insert(0, "['") + "']";\
//             }\
//             /* Static members that contain enum/key information, e.g., to map enum to string, etc. */\
//             inline static EnumMapType msEnumNames = queso::key::detail::CreateEnumMap(#KeyNames);\
//             inline static ReverseEnumMapType msEnumValues = queso::key::detail::CreateReverseEnumMap(msEnumNames);\
//             inline static std::size_t msSize = msEnumNames.size();\
//         };\
//         inline const std::string& KeyToString(EnumType value) {\
//             auto it = KesSetName##KeyType##KeyInfo::msEnumNames.find(static_cast<std::size_t>(value));\
//             if (it != KesSetName##KeyType##KeyInfo::msEnumNames.end()) {\
//                 return it->second;\
//             }\
//             QuESo_ERROR << "Invalid enum value. Possible values are: " + KesSetName##KeyType##KeyInfo::StaticGetAllKeyNames();\
//         }\
//         template<typename TType,\
//                  typename = std::enable_if_t<std::is_same<TType, KesSetName##KeyType>::value>>\
//         inline EnumType StringToKey(const std::string& rName) {\
//             const auto it = KesSetName##KeyType##KeyInfo::msEnumValues.find(rName);\
//             if (it != KesSetName##KeyType##KeyInfo::msEnumValues.end()) {\
//                 return static_cast<EnumType>(it->second);\
//             }\
//             QuESo_ERROR << "Invalid enum name. Possible values are: " + KesSetName##KeyType##KeyInfo::StaticGetAllKeyNames();\
//         }\
//         inline bool IsCorrectType(Unique<KeyInformation>& pKeyInformation, EnumType Value) {\
//             return pKeyInformation->GetKeyTypeInfo() == typeid(Value);\
//         }\
//         template<typename TType,\
//                  typename = std::enable_if_t<std::is_same<TType, EnumType>::value>>\
//         inline KesSetName##KeyType GetKeyBaseType() noexcept {\
//             return KesSetName##KeyType{};\
//         }\
//         inline std::ostream& operator<<(std::ostream& outStream, EnumType Value) {\
//             outStream << KeyToString(Value);\
//             return outStream;\
//         }\
//     }

// /**
//  * @brief Macro to register a key set for a given key set name and key type.
//  *
//  * This macro defines an enum within a struct for the given key set name and key type
//  * and creates the respective KeyInformation (QuESo_CREATE_KEY_INFO macro), which allows
//  * to get more information about the enum/key set, e.g., KeyInformation allows to map from key to string.
//  *
//  * @param KesSetName The name of the key set.
//  * @param KeyType The type of the key (possible options List, SubDict, DataSet).
//  * @param KeyNames The names of the actual keys, specified as a QuESo_LIST (QuESo_LIST macro).
//  *
//  * Example usage:
//  * @code
//  * QuESo_REGISTER_KEYS_1(KesSetName, List, QuESo_LIST(KeyName1 = 0, KeyName2, KeyName3 = 5) )
//  * @endcode
//  */
// #define QuESo_REGISTER_KEY_SET_1(KesSetName, KeyType, KeyNames) \
//     struct KesSetName {\
//         enum KeyTo##KeyType {KeyNames};\
//     };\
//     QuESo_CREATE_KEY_INFO(KesSetName, KeyType, QuESo_LIST(KeyNames), KesSetName::KeyTo##KeyType ) \

// /**
//  * @brief Macro to register a key set for a given key set name and two key types.
//  * 
//  * This macro defines two enums within a struct for the given key set name and key types
//  * and creates the respective KeyInformation (QuESo_CREATE_KEY_INFO macro), which allows
//  * to get more information about the enum/key set, e.g., KeyInformation allows to map 
//  * from key to string. 
//  * For each key type a different enum is created. This allows to distinguish between keys 
//  * that access e.g. List or DataSets. 
//  *
//  * @param KesSetName The name of the key set.
//  * @param KeyType1 The type of the first key collection (possible options List, SubDict, DataSet).
//  * @param KeyNames1 The names of the actual keys for the first key collection, specified as a QuESo_LIST (QuESo_LIST macro).
//  * @param KeyType2 The type of the second key (possible options List, SubDict, DataSet).
//  * @param KeyNames2 The names of the actual keys for the second key type, specified as a QuESo_LIST (QuESo_LIST macro).
//  *
//  * Example usage:
//  * @code
//  * QuESo_REGISTER_KEYS_2(KesSetName, List, QuESo_LIST(KeyName1 = 0, KeyName2), SubDict, QuESo_LIST(KeyName3, KeyName4 = 5) )
//  * @endcode
//  */
// #define QuESo_REGISTER_KEY_SET_2(KesSetName, KeyType1, KeyNames1, KeyType2, KeyNames2) \
//     namespace key {\
//     namespace detail {\
//         enum class CheckDuplicatedValuesOf##KesSetName {KeyNames1, KeyNames2};\
//     }\
//     }\
//     struct KesSetName {\
//         enum KeyTo##KeyType1 {KeyNames1};\
//         enum KeyTo##KeyType2 {KeyNames2};\
//     };\
//     QuESo_CREATE_KEY_INFO(KesSetName, KeyType1, QuESo_LIST(KeyNames1), KesSetName::KeyTo##KeyType1)\
//     QuESo_CREATE_KEY_INFO(KesSetName, KeyType2, QuESo_LIST(KeyNames2), KesSetName::KeyTo##KeyType2)\

// /**
//  * @brief Macro to register a key set for a given key set name and three key types.
//  * 
//  * This macro defines three enums within a struct for the given key set name and key types
//  * and creates the respective KeyInformation (QuESo_CREATE_KEY_INFO macro), which allows
//  * to get more information about the enum/key set, e.g., KeyInformation allows to map 
//  * from key to string. 
//  * For each key type a different enum is created. This allows to distinguish between keys 
//  * that access e.g. a List or a DataSet. 
//  *
//  * @param KesSetName The name of the key set.
//  * @param KeyType1 The type of the first key collection (possible options List, SubDict, DataSet).
//  * @param KeyNames1 The names of the actual keys for the first key collection, specified as a QuESo_LIST (QuESo_LIST macro).
//  * @param KeyType2 The type of the second key (possible options List, SubDict, DataSet).
//  * @param KeyNames2 The names of the actual keys for the second key type, specified as a QuESo_LIST (QuESo_LIST macro).
//  * @param KeyType3 The type of the third key (possible options List, SubDict, DataSet).
//  * @param KeyNames3 The names of the actual keys for the third key type, specified as a QuESo_LIST (QuESo_LIST macro).
//  *
//  * Example usage:
//  * @code
//  * QuESo_REGISTER_KEYS_3(KesSetName, List, QuESo_LIST(KeyName1 = 0, KeyName2), SubDict, QuESo_LIST(KeyName3, KeyName4 = 5), DataSet, QuESo_LIST(KeyName5, KeyName6 = 10) )
//  * @endcode
//  */
// #define QuESo_REGISTER_KEY_SET_3(KesSetName, KeyType1, KeyNames1, KeyType2, KeyNames2, KeyType3, KeyNames3)\
//     namespace key {\
//     namespace detail {\
//         enum class CheckDuplicatedValuesOf##KesSetName {KeyNames1, KeyNames2, KeyNames3};\
//     }\
//     }\
//     struct KesSetName {\
//         enum KeyTo##KeyType1 {KeyNames1};\
//         enum KeyTo##KeyType2 {KeyNames2};\
//         enum KeyTo##KeyType3 {KeyNames3};\
//     };\
//     QuESo_CREATE_KEY_INFO(KesSetName, KeyType1, QuESo_LIST(KeyNames1), KesSetName::KeyTo##KeyType1)\
//     QuESo_CREATE_KEY_INFO(KesSetName, KeyType2, QuESo_LIST(KeyNames2), KesSetName::KeyTo##KeyType2)\
//     QuESo_CREATE_KEY_INFO(KesSetName, KeyType3, QuESo_LIST(KeyNames3), KesSetName::KeyTo##KeyType3)\

// #endif // End KEYS_INCLUDE_HPP