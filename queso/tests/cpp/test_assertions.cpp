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

//// External includes
#include <doctest/doctest.h>

//// STL includes
#include <array>
#include <memory>
#include <string>

//// Project includes
#include "queso/includes/assert.hpp"
#include "queso/includes/parallel_exception.hpp"

namespace queso::testing {
namespace {

    constexpr bool ValidConstantExpression()
    {
        QuESo_CONTRACT_ASSERT(true, "Valid contract at compile time.");
        QuESo_DEBUG_ASSERT(true, "Valid debug assertion at compile time.");
        return true;
    }

    static_assert(ValidConstantExpression());

}  // namespace

TEST_SUITE("AssertionTestSuite")
{

    TEST_CASE("AssertionTestSuite - Contracts Always Throw With Caller Context")
    {
        const auto expected_line = __LINE__ + 2;
        try {
            QuESo_CONTRACT_ASSERT(false, "Required invariant.");
            FAIL("Contract did not throw.");
        } catch (const Exception& r_error) {
            const std::string message = r_error.what();
            CHECK(message.find("false") != std::string::npos);
            CHECK(message.find("Required invariant.") != std::string::npos);
            CHECK(message.find("test_assertions.cpp") != std::string::npos);
            CHECK(message.find("Line: " + std::to_string(expected_line)) != std::string::npos);
            CHECK(message.find(__func__) != std::string::npos);
        }
    }

    TEST_CASE("AssertionTestSuite - Conditions Evaluate Once And Messages Only On Failure")
    {
        int condition_calls = 0;
        int message_calls = 0;
        const auto condition = [&] {
            ++condition_calls;
            return true;
        };
        const auto message = [&] {
            ++message_calls;
            return "Failure message.";
        };
        // Counters intentionally observe macro evaluation semantics, not application state.
        QuESo_CONTRACT_ASSERT(condition(), message());
        CHECK((condition_calls) == (1));
        CHECK((message_calls) == (0));
        CHECK_THROWS_AS(QuESo_CONTRACT_ASSERT(!condition(), message()), Exception);
        CHECK((condition_calls) == (2));
        CHECK((message_calls) == (1));

        QuESo_DEBUG_ASSERT(condition(), message());
        CHECK((condition_calls) == (NOTDEBUG ? 2 : 3));
        CHECK((message_calls) == (1));
        const int before_diagnostic = condition_calls;
        QuESo_DIAGNOSTIC_ASSERT(condition(), message());
        CHECK((condition_calls) == (before_diagnostic + static_cast<int>(QuESo_DIAGNOSTICS_ENABLED())));
        CHECK((message_calls) == (1));
    }

    TEST_CASE("AssertionTestSuite - Disabled Checks Do Not Evaluate Failure Arguments")
    {
        int message_calls = 0;
        const auto message = [&] {
            ++message_calls;
            return "Failure message.";
        };
        if constexpr (NOTDEBUG) {
            QuESo_DEBUG_ASSERT(false, message());
        } else {
            CHECK_THROWS_AS(QuESo_DEBUG_ASSERT(false, message()), Exception);
        }
        CHECK((message_calls) == (NOTDEBUG ? 0 : 1));
        const int before_diagnostic = message_calls;
        if (QuESo_DIAGNOSTICS_ENABLED()) {
            CHECK_THROWS_AS(QuESo_DIAGNOSTIC_ASSERT(false, message()), Exception);
        } else {
            QuESo_DIAGNOSTIC_ASSERT(false, message());
        }
        CHECK((message_calls) == (before_diagnostic + static_cast<int>(QuESo_DIAGNOSTICS_ENABLED())));
    }

    TEST_CASE("AssertionTestSuite - Assertions Accept Explicit Bool And Are Single Statements")
    {
        const auto value = std::make_unique<int>(1);
        QuESo_CONTRACT_ASSERT(value, "Noncopyable explicit Boolean.");
        QuESo_DEBUG_ASSERT(value, "Noncopyable explicit Boolean.");
        QuESo_DIAGNOSTIC_ASSERT(value, "Noncopyable explicit Boolean.");
        if (value)
            QuESo_CONTRACT_ASSERT(true, "Single statement.");
        else
            FAIL("Unexpected else branch.");
        if (value)
            QuESo_DEBUG_ASSERT(true, "Single statement.");
        else
            FAIL("Unexpected else branch.");
        if (value)
            QuESo_DIAGNOSTIC_ASSERT(true, "Single statement.");
        else
            FAIL("Unexpected else branch.");
    }

    TEST_CASE("AssertionTestSuite - Worker Failure Is Rethrown After Region")
    {
        detail::ParallelException failure;
        std::array<int, 8> completed{};
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int index = 0; index < 8; ++index) {
            failure.Execute([&] {
                completed[static_cast<std::size_t>(index)] = 1;
                QuESo_CONTRACT_ASSERT(index != 3, "Worker contract failure.");
            });
        }
        for (int value : completed) { CHECK((value) == (1)); }
        CHECK_THROWS_AS(failure.Rethrow(), Exception);
        detail::ParallelException success;
        success.Execute([] {});
        CHECK_NOTHROW(success.Rethrow());
    }

}  // TEST_SUITE

}  // namespace queso::testing
