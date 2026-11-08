#pragma once

#include <string>

namespace queso::testing {

/// Initializes the data directory once, before the test runner starts. Performs no file I/O.
void InitializeTestDataDirectory(std::string Directory);

/// Returns the data directory, which remains immutable while tests execute.
[[nodiscard]] const std::string& TestDataDirectory() noexcept;

}  // namespace queso::testing
