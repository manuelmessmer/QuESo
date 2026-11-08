#include "thingi10k_options.hpp"

#include <charconv>
#include <optional>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace queso::testing {
namespace {
    std::optional<Thingi10KOptions> options;

    [[nodiscard]] IndexType ParseCount(std::string_view Value)
    {
        IndexType count{};
        const auto result = std::from_chars(Value.data(), Value.data() + Value.size(), count);
        if (result.ec != std::errc{} || result.ptr != Value.data() + Value.size() || count == 0) {
            throw std::invalid_argument("Thingi10K element counts must be positive integers.");
        }
        return count;
    }
}  // namespace

void InitializeThingi10KOptions(std::span<const std::string_view> Arguments)
{
    if (options) { throw std::logic_error("Thingi10K options are already initialized."); }
    if (Arguments.size() < 4 || Arguments.size() > 5) {
        throw std::invalid_argument("Expected: -- single/small_set/large_set input n_min n_max [model_list_directory]");
    }
    Thingi10KMode mode;
    if (Arguments[0] == "single") {
        mode = Thingi10KMode::Single;
    } else if (Arguments[0] == "small_set") {
        mode = Thingi10KMode::SmallSet;
    } else if (Arguments[0] == "large_set") {
        mode = Thingi10KMode::LargeSet;
    } else {
        throw std::invalid_argument("Unknown Thingi10K corpus mode.");
    }
    const IndexType minimum = ParseCount(Arguments[2]);
    const IndexType maximum = ParseCount(Arguments[3]);
    if (maximum < minimum) { throw std::invalid_argument("Thingi10K n_max must be at least n_min."); }
    if (Arguments[1].empty() || (mode != Thingi10KMode::Single && (Arguments.size() != 5 || Arguments[4].empty()))) {
        throw std::invalid_argument("Thingi10K input and corpus model-list directory must not be empty.");
    }
    options.emplace(
        Thingi10KOptions{ mode,
                          std::string(Arguments[1]),
                          minimum,
                          maximum,
                          Arguments.size() == 5 ? std::string(Arguments[4]) : std::string{} }
    );
}

const Thingi10KOptions& GetThingi10KOptions()
{
    if (!options) { throw std::invalid_argument("Thingi10K execution requires corpus arguments after --."); }
    return *options;
}

}  // namespace queso::testing
