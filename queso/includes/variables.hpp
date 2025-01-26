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
    template<typename TMap, typename TVector>
    inline bool CheckIfAllKeysAreSet( const TMap& rMap, const TVector& rVector ) {
        if( rMap.size() != rVector.size() ) {
            return false;
        }
        for( IndexType i = 0; i < rVector.size(); ++i ) {
            if( rMap.find(i) == rMap.end() ) {
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
    template<typename TType>\
    inline bool IsCorrectType(KeySet::KeyToDataSet Key) {\
        static const std::unordered_map<IndexType, std::type_index> type_index_map = { __VA_ARGS__ }; \
        QuESo_ASSERT( variable::detail::CheckIfAllKeysAreSet(type_index_map, key::KeySet##DataSet##KeyInfo::msEnumNames), "Types of Keys are not set correctly." );\
        const auto map_it = type_index_map.find(static_cast<std::size_t>(Key));\
        if (map_it != type_index_map.end()) {\
            return map_it->second == std::type_index(typeid(TType));\
        }\
        QuESo_ERROR << "This should be unreachable.";\
    }\
    }\

#endif // End VARIABLES_INCLUDE_HPP