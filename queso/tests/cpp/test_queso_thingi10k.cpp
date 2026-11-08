#define DOCTEST_CONFIG_IMPLEMENT
#include "helpers/thingi10k_options.hpp"
#include <doctest/doctest.h>
#include <exception>
#include <iostream>
#include <vector>

int main(int ArgumentCount, char** pArguments)
{
    for (int index = 1; index < ArgumentCount; ++index) {
        if (std::string(pArguments[index]) == "--") {
            std::vector<std::string_view> arguments;
            for (++index; index < ArgumentCount; ++index) { arguments.emplace_back(pArguments[index]); }
            try {
                if (!arguments.empty()) { queso::testing::InitializeThingi10KOptions(arguments); }
            } catch (const std::exception& r_error) {
                std::cerr << r_error.what() << '\n';
                return 2;
            }
            break;
        }
    }
    doctest::Context context(ArgumentCount, pArguments);
    return context.run();
}
