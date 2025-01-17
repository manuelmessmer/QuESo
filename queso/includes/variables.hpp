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

#ifndef VARIABLES_INCLUDE_HPP
#define VARIABLES_INCLUDE_HPP

#include<unordered_map>
#include<typeindex>

#define QuESo_VARIABLE(Key, ValueType) {Key, std::type_index(typeid(ValueType))}

namespace queso {
namespace variable {
namespace detail {
    template<typename TMap1, typename TMap2>
    inline bool CheckIfAllKeysAreSet( const TMap1& rEnumMap, const TMap2& rIndexMap ) {
        if( rEnumMap.size() != rIndexMap.size() ) {
            return false;
        }
        for (const auto& pair : rEnumMap) {
            if (rIndexMap.find(pair.first) == rIndexMap.end()) {
                return false;
            }
        }
        for (const auto& pair : rIndexMap) {
            if (rEnumMap.find(pair.first) == rEnumMap.end()) {
                return false;
            }
        }
        return true;
    }
} // End namespace detail
} // End namespace variable
} // End namesapce queso

#define QuESo_REGISTER_DATASET_VARIABLES(KeySet, ...) \
    namespace variable {\
    inline std::type_index GetValueTypeIndex(KeySet::KeyToDataSet Key) noexcept(NOTDEBUG) {\
        static const std::unordered_map<IndexType, std::type_index> type_index_map = { __VA_ARGS__ }; \
        QuESo_ASSERT( variable::detail::CheckIfAllKeysAreSet(type_index_map, key::KeySet##DataSet##KeyInfo::msEnumNames), "Types of Keys are not set correctly." );\
        const auto map_it = type_index_map.find(static_cast<std::size_t>(Key));\
        if (map_it != type_index_map.end()) {\
            return map_it->second;\
        }\
        QuESo_ERROR << "This should be unreachable.";\
    }\
    }\

#endif // End VARIABLES_INCLUDE_HPP