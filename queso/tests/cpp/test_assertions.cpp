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
#include <boost/test/unit_test.hpp>

//// STL includes
#include <array>
#include <memory>
#include <string>

//// Project includes
#include "queso/includes/assert.hpp"
#include "queso/includes/parallel_exception.hpp"

namespace queso::Testing {
namespace {

    constexpr bool ValidConstantExpression()
    {
        QuESo_CONTRACT_ASSERT(true, "Valid contract at compile time.");
        QuESo_DEBUG_ASSERT(true, "Valid debug assertion at compile time.");
        return true;
    }

    static_assert(ValidConstantExpression());

}  // namespace

BOOST_AUTO_TEST_SUITE(AssertionTestSuite)

BOOST_AUTO_TEST_CASE(ContractsAlwaysThrowWithCallerContext)
{
    const auto expected_line = __LINE__ + 2;
    try {
        QuESo_CONTRACT_ASSERT(false, "Required invariant.");
        BOOST_FAIL("Contract did not throw.");
    } catch (const Exception& r_error) {
        const std::string message = r_error.what();
        BOOST_CHECK(message.find("false") != std::string::npos);
        BOOST_CHECK(message.find("Required invariant.") != std::string::npos);
        BOOST_CHECK(message.find("test_assertions.cpp") != std::string::npos);
        BOOST_CHECK(message.find("Line: " + std::to_string(expected_line)) != std::string::npos);
        BOOST_CHECK(message.find("ContractsAlwaysThrowWithCallerContext") != std::string::npos);
    }
}

BOOST_AUTO_TEST_CASE(ConditionsEvaluateOnceAndMessagesOnlyOnFailure)
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
    BOOST_CHECK_EQUAL(condition_calls, 1);
    BOOST_CHECK_EQUAL(message_calls, 0);
    BOOST_CHECK_THROW(QuESo_CONTRACT_ASSERT(!condition(), message()), Exception);
    BOOST_CHECK_EQUAL(condition_calls, 2);
    BOOST_CHECK_EQUAL(message_calls, 1);

    QuESo_DEBUG_ASSERT(condition(), message());
    BOOST_CHECK_EQUAL(condition_calls, NOTDEBUG ? 2 : 3);
    BOOST_CHECK_EQUAL(message_calls, 1);
    const int before_diagnostic = condition_calls;
    QuESo_DIAGNOSTIC_ASSERT(condition(), message());
    BOOST_CHECK_EQUAL(condition_calls, before_diagnostic + static_cast<int>(QuESo_DIAGNOSTICS_ENABLED()));
    BOOST_CHECK_EQUAL(message_calls, 1);
}

BOOST_AUTO_TEST_CASE(DisabledChecksDoNotEvaluateFailureArguments)
{
    int message_calls = 0;
    const auto message = [&] {
        ++message_calls;
        return "Failure message.";
    };
    if constexpr (NOTDEBUG) {
        QuESo_DEBUG_ASSERT(false, message());
    } else {
        BOOST_CHECK_THROW(QuESo_DEBUG_ASSERT(false, message()), Exception);
    }
    BOOST_CHECK_EQUAL(message_calls, NOTDEBUG ? 0 : 1);
    const int before_diagnostic = message_calls;
    if (QuESo_DIAGNOSTICS_ENABLED()) {
        BOOST_CHECK_THROW(QuESo_DIAGNOSTIC_ASSERT(false, message()), Exception);
    } else {
        QuESo_DIAGNOSTIC_ASSERT(false, message());
    }
    BOOST_CHECK_EQUAL(message_calls, before_diagnostic + static_cast<int>(QuESo_DIAGNOSTICS_ENABLED()));
}

BOOST_AUTO_TEST_CASE(AssertionsAcceptExplicitBoolAndAreSingleStatements)
{
    const auto value = std::make_unique<int>(1);
    QuESo_CONTRACT_ASSERT(value, "Noncopyable explicit Boolean.");
    QuESo_DEBUG_ASSERT(value, "Noncopyable explicit Boolean.");
    QuESo_DIAGNOSTIC_ASSERT(value, "Noncopyable explicit Boolean.");
    if (value)
        QuESo_CONTRACT_ASSERT(true, "Single statement.");
    else
        BOOST_FAIL("Unexpected else branch.");
    if (value)
        QuESo_DEBUG_ASSERT(true, "Single statement.");
    else
        BOOST_FAIL("Unexpected else branch.");
    if (value)
        QuESo_DIAGNOSTIC_ASSERT(true, "Single statement.");
    else
        BOOST_FAIL("Unexpected else branch.");
}

BOOST_AUTO_TEST_CASE(WorkerFailureIsRethrownAfterRegion)
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
    for (int value : completed) { BOOST_CHECK_EQUAL(value, 1); }
    BOOST_CHECK_THROW(failure.Rethrow(), Exception);
    detail::ParallelException success;
    success.Execute([] {});
    BOOST_CHECK_NO_THROW(success.Rethrow());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
