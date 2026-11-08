#pragma once

#include <span>
#include <string>
#include <string_view>

#include "queso/includes/define.hpp"

namespace queso::testing {

/// Selects a single mesh or one of the explicit robustness corpora.
enum class Thingi10KMode { Single, SmallSet, LargeSet };

/// Typed corpus configuration parsed before test execution.
struct Thingi10KOptions
{
    Thingi10KMode mode;
    std::string input;
    IndexType minimum_elements;
    IndexType maximum_elements;
    std::string model_list_directory;
};

/// Parses corpus arguments once; does not access files or load meshes.
void InitializeThingi10KOptions(std::span<const std::string_view> Arguments);

/// Returns immutable corpus options; throws if actual corpus execution has no configuration.
[[nodiscard]] const Thingi10KOptions& GetThingi10KOptions();

}  // namespace queso::testing
