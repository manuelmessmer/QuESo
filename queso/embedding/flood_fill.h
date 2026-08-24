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

#pragma once

//// STL includes
#include <optional>
#include <vector>

//// Project includes
#include "queso/includes/define.hpp"

namespace queso::embedding {

/// @brief Cell classifications ordered by zero-based background-grid cell index.
using ElementStates = std::vector<IntersectionState>;

/// @brief Configures flood-fill execution.
struct FloodFillOptions
{
    /// @brief Number of slabs along the largest grid direction. When unset, a value is selected automatically.
    std::optional<IndexType> partition_count{};

    friend bool operator==(const FloodFillOptions&, const FloodFillOptions&) = default;
};

namespace detail {

    class DomainTopology;

    /// @brief Classifies prepared domain topology as inside, outside, or trimmed.
    /// @details Flood fill performs no geometry queries or clipping. Non-trimmed cells form components through
    ///          traversable canonical faces and are classified from exterior connectivity and directed blocked-face
    ///          votes.
    /// @param rTopology Prepared cell states, canonical face barriers, directed votes, and authoritative grid.
    /// @param Options Flood-fill partition options.
    /// @return States ordered by zero-based grid cell index.
    [[nodiscard]] ElementStates FloodFill(const DomainTopology& rTopology, FloodFillOptions Options = {});

}  // namespace detail

}  // namespace queso::embedding
