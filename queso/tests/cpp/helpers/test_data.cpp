#include "test_data.hpp"

#include <stdexcept>
#include <utility>

namespace queso::testing {
namespace {
    std::string data_directory = QUESO_TEST_DATA_DIR;
    bool is_initialized = false;
}  // namespace

void InitializeTestDataDirectory(std::string Directory)
{
    if (is_initialized) { throw std::logic_error("Test data directory is already initialized."); }
    if (Directory.empty()) { throw std::invalid_argument("Test data directory must not be empty."); }
    data_directory = std::move(Directory);
    is_initialized = true;
}

const std::string& TestDataDirectory() noexcept
{ return data_directory; }

}  // namespace queso::testing
