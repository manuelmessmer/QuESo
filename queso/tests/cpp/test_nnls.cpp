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

//// STL includes
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

//// External includes
#include <boost/test/unit_test.hpp>

//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/includes/logger.hpp"
#include "queso/solvers/nnls.h"
#include "queso/utilities/math_utilities.hpp"

namespace queso::Testing {
namespace {

    long double ResidualNorm(std::span<const double> Lhs, std::span<const double> Rhs, std::span<const double> Solution)
    {
        long double norm = 0.0L;
        for (IndexType i = 0; i < Rhs.size(); ++i) {
            long double residual = Rhs[i];
            for (IndexType j = 0; j < Solution.size(); ++j) {
                residual -= static_cast<long double>(Lhs[j * Rhs.size() + i]) * Solution[j];
            }
            norm = std::hypot(norm, residual);
        }
        return norm;
    }

    void CheckOptimality(
        std::span<const double> Lhs,
        std::span<const double> Rhs,
        std::span<const double> Solution,
        const nnls::Result& rResult
    )
    {
        BOOST_REQUIRE(rResult.status == nnls::Status::Converged);
        std::vector<long double> residual(Rhs.begin(), Rhs.end());
        long double norm = 0.0L;
        long double rhs_norm = 0.0L;
        for (IndexType i = 0; i < Rhs.size(); ++i) {
            rhs_norm = std::hypot(rhs_norm, static_cast<long double>(Rhs[i]));
            for (IndexType j = 0; j < Solution.size(); ++j) {
                residual[i] -= static_cast<long double>(Lhs[j * Rhs.size() + i]) * Solution[j];
            }
            norm = std::hypot(norm, residual[i]);
        }
        if (rhs_norm == 0.0L) {
            BOOST_CHECK_EQUAL(norm, 0.0L);
            BOOST_CHECK_EQUAL(rResult.residual_norm, 0.0);
        } else {
            BOOST_CHECK_SMALL((static_cast<long double>(rResult.residual_norm) - norm) / rhs_norm, 2e-12L);
            for (auto& r_value : residual) { r_value /= rhs_norm; }
        }
        IndexType support = 0;
        for (IndexType j = 0; j < Solution.size(); ++j) {
            BOOST_REQUIRE(std::isfinite(Solution[j]));
            BOOST_CHECK_GE(Solution[j], 0.0);
            support += Solution[j] > 0.0 ? 1 : 0;
            const auto column = Lhs.subspan(j * Rhs.size(), Rhs.size());
            long double column_scale = 0.0L;
            for (const double value : column) {
                column_scale = std::max(column_scale, std::abs(static_cast<long double>(value)));
            }
            if (column_scale == 0.0L) { continue; }
            long double dual = 0.0L;
            long double column_norm = 0.0L;
            for (IndexType i = 0; i < Rhs.size(); ++i) {
                const long double value = column[i] / column_scale;
                dual += value * residual[i];
                column_norm = std::hypot(column_norm, value);
            }
            dual /= column_norm;
            BOOST_CHECK_LE(dual, 2e-11L);
            if (Solution[j] > 0.0) { BOOST_CHECK_SMALL(dual, 2e-11L); }
        }
        BOOST_CHECK_LE(support, std::min(Rhs.size(), Solution.size()));
    }

    // Independent tiny-system oracle: enumerate all supports and solve their normal equations in long double.
    // This is deliberately test-only, and used for small, moderately scaled random integer matrices.
    long double ExhaustiveResidual(std::span<const double> Lhs, std::span<const double> Rhs, IndexType Columns)
    {
        long double best = ResidualNorm(Lhs, Rhs, std::vector<double>(Columns, 0.0));
        for (IndexType mask = 1; mask < (IndexType{ 1 } << Columns); ++mask) {
            std::vector<IndexType> active;
            for (IndexType j = 0; j < Columns; ++j) {
                if ((mask & (IndexType{ 1 } << j)) != 0) { active.push_back(j); }
            }
            const IndexType count = active.size();
            if (count > Rhs.size()) { continue; }
            std::vector<long double> system(count * (count + 1), 0.0L);
            for (IndexType i = 0; i < count; ++i) {
                for (IndexType j = 0; j <= count; ++j) {
                    for (IndexType k = 0; k < Rhs.size(); ++k) {
                        system[i * (count + 1) + j] += static_cast<long double>(Lhs[active[i] * Rhs.size() + k])
                                                       * (j == count ? Rhs[k] : Lhs[active[j] * Rhs.size() + k]);
                    }
                }
            }
            bool is_singular = false;
            for (IndexType pivot = 0; pivot < count; ++pivot) {
                IndexType selected = pivot;
                for (IndexType i = pivot + 1; i < count; ++i) {
                    if (std::abs(system[i * (count + 1) + pivot]) > std::abs(system[selected * (count + 1) + pivot])) {
                        selected = i;
                    }
                }
                if (std::abs(system[selected * (count + 1) + pivot]) < 1e-15L) {
                    is_singular = true;
                    break;
                }
                for (IndexType j = pivot; j <= count; ++j) {
                    std::swap(system[pivot * (count + 1) + j], system[selected * (count + 1) + j]);
                }
                for (IndexType i = pivot + 1; i < count; ++i) {
                    const long double factor = system[i * (count + 1) + pivot] / system[pivot * (count + 1) + pivot];
                    for (IndexType j = pivot; j <= count; ++j) {
                        system[i * (count + 1) + j] -= factor * system[pivot * (count + 1) + j];
                    }
                }
            }
            if (is_singular) { continue; }
            std::vector<double> solution(Columns, 0.0);
            bool is_feasible = true;
            for (IndexType remaining = count; remaining != 0; --remaining) {
                const IndexType i = remaining - 1;
                long double value = system[i * (count + 1) + count];
                for (IndexType j = i + 1; j < count; ++j) {
                    value -= system[i * (count + 1) + j] * solution[active[j]];
                }
                value /= system[i * (count + 1) + i];
                if (value < -1e-14L) { is_feasible = false; }
                solution[active[i]] = static_cast<double>(std::max(0.0L, value));
            }
            if (is_feasible) { best = std::min(best, ResidualNorm(Lhs, Rhs, solution)); }
        }
        return best;
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(NNLSTestSuite)

BOOST_AUTO_TEST_CASE(NNLS_Rectangle_Matrix)
{
    QuESo_INFO << "Testing :: Test nnls :: Rectangle Matrix" << std::endl;

    std::vector<double> a = { 0.1, 0.3, 0.2, 0.5, 0.2, 0.7, 0.5, 0.3, 3.0, 0.6, 1.1, 0.9 };

    std::vector<double> rhs = { 0.1, 0.1, 0.7, 0.3 };

    std::vector<double> solution(3);
    nnls::Workspace workspace;
    const auto result = nnls::Solve({ .lhs = a, .rhs = rhs, .solution = solution, .r_workspace = workspace });
    BOOST_REQUIRE(result.status == nnls::Status::Converged);
    const double rnorm = result.residual_norm;

    std::vector<double> reference(3);
    reference[0] = 0.26527761625544183;
    reference[1] = 0.39411741902728925;
    reference[2] = 0.032491624806329486;

    QuESo_CHECK_NEAR(rnorm, 0.5080235032184826, 1e-10);
    for (IndexType i = 0; i < reference.size(); ++i) QuESo_CHECK_NEAR(solution[i], reference[i], 1e-10);
}

BOOST_AUTO_TEST_CASE(NNLS_Square_Matrix)
{
    QuESo_INFO << "Testing :: Test nnls :: Square Matrix" << std::endl;

    std::vector<double> a = { 0.2, 0.3, 0.87, 0.22, 0.5, 0.7, 0.5, 0.45, 6.0, 0.6, 1.1, 0.9, 0.3, 0.3, 0.2, 1.1 };

    std::vector<double> rhs = { 0.33, 0.12, 0.12, 0.3 };

    std::vector<double> solution(4);
    nnls::Workspace workspace;
    const auto result = nnls::Solve({ .lhs = a, .rhs = rhs, .solution = solution, .r_workspace = workspace });
    BOOST_REQUIRE(result.status == nnls::Status::Converged);
    const double rnorm = result.residual_norm;

    std::vector<double> reference(4);
    reference[0] = 0.015692959503716887;
    reference[1] = 0.0347636417196625;
    reference[2] = 0.0404670415359823;
    reference[3] = 0.22225779341177274;

    QuESo_CHECK_NEAR(rnorm, 0.0, 1e-10);
    for (IndexType i = 0; i < reference.size(); ++i) QuESo_CHECK_NEAR(solution[i], reference[i], 1e-10);
}

BOOST_AUTO_TEST_CASE(BoundaryAndRankDeficientSystems)
{
    nnls::Workspace workspace;
    const std::vector<std::vector<double>> matrices{
        { 1.0, 0.0, 0.0, 1.0 },
        { 2.0, 1.0, 1.0, 1.0 },  // First variable must leave the active set.
        { 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0 },
    };
    const std::vector<double> rhs{ -0.1, 1.0 };
    for (const auto& r_matrix : matrices) {
        auto lhs = r_matrix;
        auto rhs_buffer = rhs;
        std::vector<double> solution(lhs.size() / rhs.size(), -7.0);
        const auto result =
            nnls::Solve({ .lhs = lhs, .rhs = rhs_buffer, .solution = solution, .r_workspace = workspace });
        CheckOptimality(r_matrix, rhs, solution, result);
    }
}

BOOST_AUTO_TEST_CASE(RandomSystemsAgreeWithExhaustiveOracle)
{
    nnls::Workspace workspace;
    std::uint32_t state = 71023;
    const auto next = [&state]() {
        state = state * 1664525U + 1013904223U;
        return static_cast<double>((state >> 16U) % 11U) - 5.0;
    };
    for (IndexType case_index = 0; case_index < 160; ++case_index) {
        const IndexType rows = 1 + case_index % 5;
        const IndexType columns = 1 + case_index % 7;
        std::vector<double> lhs(rows * columns);
        std::vector<double> rhs(rows);
        std::generate(lhs.begin(), lhs.end(), next);
        std::generate(rhs.begin(), rhs.end(), next);
        const auto original_lhs = lhs;
        const auto original_rhs = rhs;
        std::vector<double> solution(columns);
        const auto result = nnls::Solve({ .lhs = lhs, .rhs = rhs, .solution = solution, .r_workspace = workspace });
        BOOST_TEST_CONTEXT("case " << case_index)
        {
            CheckOptimality(original_lhs, original_rhs, solution, result);
            BOOST_CHECK_SMALL(
                ResidualNorm(original_lhs, original_rhs, solution)
                    - ExhaustiveResidual(original_lhs, original_rhs, columns),
                3e-11L
            );
        }
    }
}

BOOST_AUTO_TEST_CASE(ScaledColumnsAndPhysicalWeights)
{
    nnls::Workspace workspace;
    for (const double scale : { 1e-100, 1.0, 1e100 }) {
        const std::vector<double> original_lhs{ 1.0, 0.0, 0.0, 1e-20 };
        const std::vector<double> original_rhs{ 1e-20 * scale, scale };
        auto lhs = original_lhs;
        auto rhs = original_rhs;
        std::vector<double> solution(2);
        const auto result = nnls::Solve({ .lhs = lhs, .rhs = rhs, .solution = solution, .r_workspace = workspace });
        CheckOptimality(original_lhs, original_rhs, solution, result);
        BOOST_CHECK_CLOSE_FRACTION(solution[0], original_rhs[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(solution[1], scale * 1e20, 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(IterationLimitReturnsActualFeasibleResidual)
{
    const std::vector<double> original_lhs{ 2.0, 1.0, 1.0, 1.0 };
    const std::vector<double> original_rhs{ 0.0, 1.0 };
    auto lhs = original_lhs;
    auto rhs = original_rhs;
    std::vector<double> solution(2);
    nnls::Workspace workspace;
    const auto result = nnls::Solve(
        { .lhs = lhs, .rhs = rhs, .solution = solution, .r_workspace = workspace, .options = { .max_iterations = 1 } }
    );
    BOOST_CHECK(result.status == nnls::Status::IterationLimit);
    BOOST_CHECK_EQUAL(result.iterations, 1);
    BOOST_CHECK_GT(result.residual_norm, 0.8);
    BOOST_CHECK_SMALL(
        static_cast<long double>(result.residual_norm) - ResidualNorm(original_lhs, original_rhs, solution), 1e-14L
    );
    for (const double value : solution) { BOOST_CHECK_GE(value, 0.0); }
}

BOOST_AUTO_TEST_CASE(RankRejectionIsNotConvergence)
{
    std::vector<double> lhs{ 2.0, 0.0, 1.0, 1e-3 };
    std::vector<double> rhs{ 1.0, 1.0 };
    std::vector<double> solution(2);
    nnls::Workspace workspace;
    const auto result = nnls::Solve(
        { .lhs = lhs, .rhs = rhs, .solution = solution, .r_workspace = workspace, .options = { .rank_tolerance = 0.1 } }
    );
    BOOST_CHECK(result.status == nnls::Status::NumericalFailure);
    BOOST_CHECK_SMALL(result.residual_norm - 1.0, 1e-14);
}

BOOST_AUTO_TEST_CASE(StableNormsAndArithmeticFailure)
{
    for (const double scale : { 1e-200, 1e200 }) {
        const std::vector<double> values{ 3.0 * scale, 4.0 * scale };
        BOOST_CHECK_CLOSE_FRACTION(Math::L2Norm(values), 5.0 * scale, 1e-15);
    }
    std::vector<double> lhs{ 1e-200 };
    std::vector<double> rhs{ 1e-100 };
    std::vector<double> solution(1);
    nnls::Workspace workspace;
    auto result = nnls::Solve({ .lhs = lhs, .rhs = rhs, .solution = solution, .r_workspace = workspace });
    BOOST_REQUIRE(result.status == nnls::Status::Converged);
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1e100, 1e-15);

    lhs[0] = 1.0;
    rhs[0] = std::numeric_limits<double>::max();
    result = nnls::Solve({ .lhs = lhs, .rhs = rhs, .solution = solution, .r_workspace = workspace });
    BOOST_CHECK(result.status == nnls::Status::NumericalFailure);
    BOOST_CHECK(std::isinf(result.residual_norm));  // Explicitly unavailable, never a fabricated zero.
}

BOOST_AUTO_TEST_CASE(NormalizationRangeBoundaries)
{
    const double tiny = std::numeric_limits<double>::denorm_min();
    constexpr double minimum = std::numeric_limits<double>::min();
    struct RangeCase
    {
        std::vector<double> lhs;
        std::vector<double> rhs;
        nnls::Status status;
        nnls::Options options{};
    };
    const std::vector<RangeCase> cases{
        { { tiny, tiny }, { 1e-300, 0.0 }, nnls::Status::NumericalFailure },
        { { 1.0, 0.0, 0.0, 1.0 }, { tiny, tiny }, nnls::Status::NumericalFailure },
        { { 1.0, 0.0, 0.0, 1.0 }, { 1.0, tiny }, nnls::Status::NumericalFailure },
        { { 1.0, 0.0, 1.0, 2.0 * tiny },
          { 1.0, 1e-300 },
          nnls::Status::NumericalFailure,
          { .dual_tolerance = tiny, .rank_tolerance = tiny } },
        { { std::nextafter(minimum, 0.0) }, { 1.0 }, nnls::Status::NumericalFailure },
        { { minimum, 0.0, 0.0, minimum }, { minimum, minimum }, nnls::Status::Converged },
        { { std::nextafter(minimum, 1.0) }, { minimum }, nnls::Status::Converged },
        { { tiny, tiny }, { 0.0, 0.0 }, nnls::Status::Converged },
        { {}, { tiny }, nnls::Status::Converged },
    };
    nnls::Workspace workspace;
    for (IndexType i = 0; i < cases.size(); ++i) {
        BOOST_TEST_CONTEXT("range case " << i)
        {
            auto lhs = cases[i].lhs;
            auto rhs = cases[i].rhs;
            std::vector<double> solution(lhs.size() / rhs.size());
            const auto result = nnls::Solve(
                { .lhs = lhs, .rhs = rhs, .solution = solution, .r_workspace = workspace, .options = cases[i].options }
            );
            BOOST_REQUIRE(result.status == cases[i].status);
            if (result.status == nnls::Status::Converged) {
                CheckOptimality(cases[i].lhs, cases[i].rhs, solution, result);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(KnownOptimaAcrossScalesAndPermutations)
{
    nnls::Workspace workspace;
    for (const IndexType rows : { IndexType{ 27 }, IndexType{ 64 } }) {
        for (const double delta : { 0.0, std::ldexp(1.0, -12), 1.0 }) {
            for (const int exponent : { -1000, 0, 400 }) {
                for (const bool is_permuted : { false, true }) {
                    BOOST_TEST_CONTEXT(
                        "rows " << rows << ", delta " << delta << ", exponent " << exponent << ", permuted "
                                << is_permuted
                    )
                    {
                        const IndexType columns = rows + 2;
                        const double scale = std::ldexp(1.0, exponent);
                        std::vector<double> lhs(rows * columns, 0.0), rhs(rows, 0.0), solution(columns);
                        const auto entry = [&](IndexType I, IndexType J) -> double& {
                            return lhs[(is_permuted ? columns - 1 - J : J) * rows + (is_permuted ? rows - 1 - I : I)];
                        };
                        // x=(1,2,0,...) has residual e_last. Every inactive column has strictly negative dual.
                        entry(0, 0) = entry(0, 1) = scale;
                        entry(1, 1) = delta * scale;
                        for (IndexType j = 2; j < columns; ++j) {
                            entry(2 + (j - 2) % (rows - 3), j) = scale;
                            entry(rows - 1, j) = -scale;
                        }
                        rhs[0] = 3.0 * scale;
                        rhs[1] = 2.0 * delta * scale;
                        rhs.back() = scale;
                        if (is_permuted) { std::ranges::reverse(rhs); }
                        const auto original_lhs = lhs;
                        const auto original_rhs = rhs;
                        const auto result =
                            nnls::Solve({ .lhs = lhs, .rhs = rhs, .solution = solution, .r_workspace = workspace });
                        CheckOptimality(original_lhs, original_rhs, solution, result);
                        BOOST_CHECK_SMALL(result.residual_norm / scale - 1.0, 2e-12);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(InvalidInputsAreNotModified)
{
    nnls::Workspace workspace;
    const std::vector<double> original_rhs{ 1.0, 2.0 };
    std::vector<double> lhs{ 1.0, 2.0, 3.0 };
    auto rhs = original_rhs;
    std::vector<double> solution(2, 7.0);
    auto result = nnls::Solve({ .lhs = lhs, .rhs = rhs, .solution = solution, .r_workspace = workspace });
    BOOST_CHECK(result.status == nnls::Status::InvalidInput);
    BOOST_CHECK_EQUAL(solution[0], 7.0);
    BOOST_CHECK(rhs == original_rhs);
    lhs = { 1.0, 0.0, 0.0, 1.0 };
    result = nnls::Solve(
        { .lhs = lhs,
          .rhs = rhs,
          .solution = solution,
          .r_workspace = workspace,
          .options = { .dual_tolerance = -1.0 } }
    );
    BOOST_CHECK(result.status == nnls::Status::InvalidInput);
    for (const double invalid : { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN() }) {
        lhs[0] = invalid;
        result = nnls::Solve({ .lhs = lhs, .rhs = rhs, .solution = solution, .r_workspace = workspace });
        BOOST_CHECK(result.status == nnls::Status::InvalidInput);
        BOOST_CHECK_EQUAL(solution[0], 7.0);
    }
    lhs[0] = 1.0;
    result =
        nnls::Solve({ .lhs = lhs, .rhs = rhs, .solution = std::span<double>(lhs).first(2), .r_workspace = workspace });
    BOOST_CHECK(result.status == nnls::Status::InvalidInput);
    result = nnls::Solve({ .lhs = lhs, .rhs = rhs, .solution = rhs, .r_workspace = workspace });
    BOOST_CHECK(result.status == nnls::Status::InvalidInput);
    result = nnls::Solve(
        { .lhs = lhs, .rhs = std::span<double>(lhs).subspan(1, 2), .solution = solution, .r_workspace = workspace }
    );
    BOOST_CHECK(result.status == nnls::Status::InvalidInput);
}

BOOST_AUTO_TEST_CASE(EmptySystemsAndZeroRhs)
{
    nnls::Workspace workspace;
    std::vector<double> empty;
    std::vector<double> values{ 3.0, 4.0 };
    auto result = nnls::Solve({ .lhs = empty, .rhs = values, .solution = empty, .r_workspace = workspace });
    BOOST_CHECK(result.status == nnls::Status::Converged);
    BOOST_CHECK_EQUAL(result.residual_norm, 5.0);
    result = nnls::Solve({ .lhs = empty, .rhs = empty, .solution = values, .r_workspace = workspace });
    BOOST_CHECK(result.status == nnls::Status::Converged);
    BOOST_CHECK_EQUAL(values[0], 0.0);
    result = nnls::Solve({ .lhs = empty, .rhs = empty, .solution = empty, .r_workspace = workspace });
    BOOST_CHECK(result.status == nnls::Status::Converged);
    std::vector<double> lhs{ 1.0, 2.0, 3.0, 4.0 };
    std::vector<double> rhs(2, 0.0);
    result = nnls::Solve({ .lhs = lhs, .rhs = rhs, .solution = values, .r_workspace = workspace });
    BOOST_CHECK(result.status == nnls::Status::Converged);
    BOOST_CHECK_EQUAL(result.residual_norm, 0.0);
    BOOST_CHECK_EQUAL(values[1], 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace queso::Testing
