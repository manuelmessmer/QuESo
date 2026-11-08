#define DOCTEST_CONFIG_IMPLEMENT
#include "helpers/test_data.hpp"
#include <doctest/doctest.h>
#include <exception>
#include <iostream>
#include <utility>

int main(int ArgumentCount, char** pArguments)
{
    std::string data_directory = QUESO_TEST_DATA_DIR;
    for (int index = 1; index < ArgumentCount; ++index) {
        const std::string argument = pArguments[index];
        if (argument.starts_with("--data-dir=")) { data_directory = argument.substr(11); }
    }
    try {
        queso::testing::InitializeTestDataDirectory(std::move(data_directory));
    } catch (const std::exception& r_error) {
        std::cerr << r_error.what() << '\n';
        return 2;
    }
    doctest::Context context(ArgumentCount, pArguments);
    return context.run();
}
