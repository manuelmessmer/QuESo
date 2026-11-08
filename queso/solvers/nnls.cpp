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
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>

//// Project includes
#include "queso/solvers/nnls.h"
#include "queso/utilities/math_utilities.hpp"

namespace queso::nnls {
namespace {

    // Independent accumulation chains expose portable SIMD without enabling reassociation via fast-math.
    // Both inputs have the same length. Callers supply already validated subspans.
    constexpr double Dot(std::span<const double> Left, std::span<const double> Right) noexcept
    {
        std::array<double, 4> sums{};
        const IndexType block_end = Left.size() - Left.size() % sums.size();
        IndexType i = 0;
        for (; i < block_end; i += sums.size()) {
            sums[0] += Left[i] * Right[i];
            sums[1] += Left[i + 1] * Right[i + 1];
            sums[2] += Left[i + 2] * Right[i + 2];
            sums[3] += Left[i + 3] * Right[i + 3];
        }
        double sum = (sums[0] + sums[1]) + (sums[2] + sums[3]);
        for (; i < Left.size(); ++i) { sum += Left[i] * Right[i]; }
        return sum;
    }

    bool Overlaps(std::span<const double> Left, std::span<const double> Right) noexcept
    {
        if (Left.empty() || Right.empty()) { return false; }
        // std::less supplies a total pointer order without subtracting unrelated pointers.
        const std::less<const double*> less;
        return less(Left.data(), Right.data() + Right.size()) && less(Right.data(), Left.data() + Left.size());
    }

    bool IsValidTolerance(double Value) noexcept
    { return std::isfinite(Value) && Value >= 0.0 && Value < 1.0; }

    bool IsUsableNorm(double Value) noexcept
    { return std::isfinite(Value) && (Value == 0.0 || Value >= std::numeric_limits<double>::min()); }

    // All matrix columns represent Q*A, not packed reflectors. Active columns have exact zero subdiagonals.
    // Only the trial vector may be infeasible. The caller's solution remains feasible throughout.
    class Solver
    {
    public:
        Solver(
            Parameters Arguments,
            std::span<double> Dual,
            std::span<double> ColumnNorms,
            std::span<double> Trial,
            std::span<double> Reflector,
            std::span<IndexType> Indices
        ) noexcept
            : mLhs(Arguments.lhs), mRhs(Arguments.rhs), mSolution(Arguments.solution), mDual(Dual),
              mColumnNorms(ColumnNorms), mTrial(Trial), mReflector(Reflector), mIndices(Indices), mRows(mRhs.size()),
              mColumns(mSolution.size())
        {
            constexpr double epsilon = std::numeric_limits<double>::epsilon();
            const double dimension = static_cast<double>(std::max(IndexType{ 1 }, mRows));
            mDualTolerance =
                Arguments.options.dual_tolerance == 0.0 ? 8.0 * epsilon * dimension : Arguments.options.dual_tolerance;
            mRankTolerance =
                Arguments.options.rank_tolerance == 0.0 ? epsilon * dimension : Arguments.options.rank_tolerance;
            constexpr IndexType max = std::numeric_limits<IndexType>::max();
            mMaxIterations = Arguments.options.max_iterations != 0 ? Arguments.options.max_iterations
                                                                   : (mColumns > max / 5 ? max : 5 * mColumns);
        }

        Result Run() noexcept
        {
            mRhsNorm = Math::L2Norm(mRhs);
            if (!std::isfinite(mRhsNorm)) { return Failure(); }
            std::ranges::fill(mSolution, 0.0);
            std::iota(mIndices.begin(), mIndices.end(), IndexType{ 0 });
            if (mRows == 0 || mColumns == 0 || mRhsNorm == 0.0) { return Finish(Status::Converged); }
            if (!IsUsableNorm(mRhsNorm)) { return Failure(); }
            for (IndexType j = 0; j < mColumns; ++j) {
                mColumnNorms[j] = Math::L2Norm(Column(j));
                if (!IsUsableNorm(mColumnNorms[j])) { return Failure(); }
            }

            while (true) {
                // Normalize by the current trailing residual, preserving pivot ordering and allowing recovery of
                // small but resolvable positive weights after larger residual components have been fitted.
                const double residual_norm = Math::L2Norm(mRhs.subspan(mActive));
                if (!IsUsableNorm(residual_norm)) { return Failure(); }
                if (residual_norm == 0.0) { return Finish(Status::Converged); }
                for (IndexType i = mActive; i < mRows; ++i) { mReflector[i] = mRhs[i] / residual_norm; }
                for (IndexType position = mActive; position < mColumns; ++position) {
                    const IndexType j = mIndices[position];
                    const double dual = Dot(Column(j).subspan(mActive), mReflector.subspan(mActive));
                    if (!std::isfinite(dual)) { return Failure(); }
                    mDual[j] = dual;
                }

                bool has_blocked_candidate = false;
                while (true) {
                    IndexType candidate = mColumns;
                    double largest_dual = 0.0;
                    bool has_dual_violation = false;
                    for (IndexType position = mActive; position < mColumns; ++position) {
                        const IndexType j = mIndices[position];
                        has_dual_violation |= mColumnNorms[j] > 0.0 && mDual[j] / mColumnNorms[j] > mDualTolerance;
                        // Preserve raw-dual pivot ordering, even if a larger raw dual has a smaller relative value.
                        if (mDual[j] > largest_dual) {
                            candidate = position;
                            largest_dual = mDual[j];
                        }
                    }
                    if (!has_dual_violation) {
                        return Finish(has_blocked_candidate ? Status::NumericalFailure : Status::Converged);
                    }
                    const IndexType j = mIndices[candidate];
                    const double tail_norm = Math::L2Norm(Column(j).subspan(mActive));
                    if (!IsUsableNorm(tail_norm)) { return Failure(); }
                    if (tail_norm / mColumnNorms[j] <= mRankTolerance) {
                        has_blocked_candidate |= mDual[j] / mColumnNorms[j] > mDualTolerance;
                        mDual[j] = 0.0;
                        continue;
                    }

                    const double diagonal = -std::copysign(tail_norm, Column(j)[mActive]);
                    const double first = Column(j)[mActive] / tail_norm;
                    const double divisor = first + std::copysign(1.0, first);
                    const double tau = 1.0 + std::abs(first);
                    mReflector[mActive] = 1.0;
                    for (IndexType i = mActive + 1; i < mRows; ++i) {
                        mReflector[i] = (Column(j)[i] / tail_norm) / divisor;
                    }
                    std::ranges::copy(mRhs, mTrial.begin());
                    if (!ApplyReflector(mTrial, tau)) { return Failure(); }
                    const double proposed = mTrial[mActive] / diagonal;
                    if (!std::isfinite(proposed)) { return Failure(); }
                    if (proposed <= 0.0) {
                        has_blocked_candidate |= mDual[j] / mColumnNorms[j] > mDualTolerance;
                        mDual[j] = 0.0;
                        continue;
                    }

                    // Commit only after independence and positivity checks. Rejected candidates leave Q*A and Q*b
                    // intact.
                    std::swap(mIndices[mActive], mIndices[candidate]);
                    for (IndexType position = mActive + 1; position < mColumns; ++position) {
                        if (!ApplyReflector(Column(mIndices[position]), tau)) { return Failure(); }
                    }
                    std::ranges::copy(mTrial, mRhs.begin());
                    Column(j)[mActive] = diagonal;
                    for (IndexType i = mActive + 1; i < mRows; ++i) { Column(j)[i] = 0.0; }
                    ++mActive;
                    break;
                }

                while (true) {
                    if (mIterations == mMaxIterations) { return Finish(Status::IterationLimit); }
                    ++mIterations;
                    if (!SolveTriangular()) { return Finish(Status::NumericalFailure); }
                    double alpha = 1.0;
                    IndexType boundary = mActive;
                    for (IndexType position = 0; position < mActive; ++position) {
                        if (mTrial[position] <= 0.0) {
                            const double current = mSolution[mIndices[position]];
                            const double negative_trial = -mTrial[position];
                            // x / (x-z), without overflowing x-z or dividing 0 by 0.
                            const double ratio =
                                current == 0.0 ? 0.0
                                               : (current > negative_trial
                                                      ? 1.0 / (1.0 + negative_trial / current)
                                                      : (current / negative_trial) / (1.0 + current / negative_trial));
                            if (boundary == mActive || ratio < alpha) {
                                alpha = ratio;
                                boundary = position;
                            }
                        }
                    }
                    if (boundary == mActive) {
                        for (IndexType position = 0; position < mActive; ++position) {
                            mSolution[mIndices[position]] = mTrial[position];
                        }
                        break;
                    }
                    for (IndexType position = 0; position < mActive; ++position) {
                        double& r_value = mSolution[mIndices[position]];
                        r_value = std::max(0.0, (1.0 - alpha) * r_value + alpha * mTrial[position]);
                    }
                    mSolution[mIndices[boundary]] = 0.0;
                    for (IndexType position = 0; position < mActive;) {
                        if (mSolution[mIndices[position]] == 0.0) {
                            if (!RemoveColumn(position)) { return Failure(); }
                        } else {
                            ++position;
                        }
                    }
                }
            }
        }

    private:
        std::span<double> Column(IndexType J) const noexcept
        { return mLhs.subspan(J * mRows, mRows); }

        bool ApplyReflector(std::span<double> Values, double Tau) const noexcept;
        bool SolveTriangular() noexcept;
        bool RemoveColumn(IndexType Position) noexcept;

        Result Finish(Status SolveStatus) noexcept
        {
            if (!std::ranges::all_of(mLhs, [](double Value) { return std::isfinite(Value); })) { return Failure(); }
            // Include the leading block: on exhaustion x need not solve the most recently expanded triangular system.
            std::ranges::copy(mRhs, mTrial.begin());
            for (IndexType position = 0; position < mActive; ++position) {
                const IndexType j = mIndices[position];
                for (IndexType i = 0; i <= position; ++i) { mTrial[i] -= Column(j)[i] * mSolution[j]; }
            }
            if (!std::ranges::all_of(mTrial, [](double Value) { return std::isfinite(Value); })) { return Failure(); }
            const double residual = Math::L2Norm(mTrial);
            if (!std::isfinite(residual)) { return Failure(); }
            return { .status = SolveStatus, .residual_norm = residual, .iterations = mIterations };
        }

        Result Failure() const noexcept
        { return { .status = Status::NumericalFailure, .iterations = mIterations }; }

        std::span<double> mLhs;
        std::span<double> mRhs;
        std::span<double> mSolution;
        std::span<double> mDual;
        std::span<double> mColumnNorms;
        std::span<double> mTrial;
        std::span<double> mReflector;
        std::span<IndexType> mIndices;
        IndexType mRows;
        IndexType mColumns;
        IndexType mActive = 0;
        IndexType mIterations = 0;
        IndexType mMaxIterations;
        double mRhsNorm = 0.0;
        double mDualTolerance;
        double mRankTolerance;
    };

    bool Solver::ApplyReflector(std::span<double> Values, double Tau) const noexcept
    {
        const auto reflector = mReflector;
        const double projection = Dot(Values.subspan(mActive), reflector.subspan(mActive));
        const double factor = Tau * projection;
        if (!std::isfinite(factor)) { return false; }
        for (IndexType i = mActive; i < mRows; ++i) { Values[i] -= factor * reflector[i]; }
        // Keep updates vectorizable. Dual/triangular work and Finish detect nonfinite entries before success.
        return true;
    }

    bool Solver::SolveTriangular() noexcept
    {
        const auto trial = mTrial;
        std::copy_n(mRhs.begin(), mActive, trial.begin());
        for (IndexType remaining = mActive; remaining != 0; --remaining) {
            const IndexType j = remaining - 1;
            const auto values = Column(mIndices[j]);
            const double coefficient = trial[j] / values[j];
            if (!std::isfinite(coefficient)) { return false; }
            trial[j] = coefficient;
            for (IndexType i = 0; i < j; ++i) { trial[i] -= values[i] * coefficient; }
        }
        return true;
    }

    bool Solver::RemoveColumn(IndexType Position) noexcept
    {
        const IndexType removed = mIndices[Position];
        for (IndexType i = Position; i + 1 < mActive; ++i) {
            const IndexType j = mIndices[i + 1];
            const double upper = Column(j)[i];
            const double lower = Column(j)[i + 1];
            const double radius = std::hypot(upper, lower);
            if (!std::isfinite(radius)) { return false; }
            const double cosine = radius == 0.0 ? 1.0 : upper / radius;
            const double sine = radius == 0.0 ? 0.0 : lower / radius;
            for (IndexType k = 0; k < mColumns; ++k) {
                auto values = Column(k);
                const double first = values[i];
                values[i] = cosine * first + sine * values[i + 1];
                values[i + 1] = -sine * first + cosine * values[i + 1];
            }
            Column(j)[i] = radius;
            Column(j)[i + 1] = 0.0;
            const double first = mRhs[i];
            mRhs[i] = cosine * first + sine * mRhs[i + 1];
            mRhs[i + 1] = -sine * first + cosine * mRhs[i + 1];
            if (!std::isfinite(mRhs[i]) || !std::isfinite(mRhs[i + 1])) { return false; }
            mIndices[i] = j;
        }
        --mActive;
        mIndices[mActive] = removed;
        return true;
    }

}  // namespace

Result Solve(Parameters Arguments)
{
    const IndexType rows = Arguments.rhs.size();
    const IndexType columns = Arguments.solution.size();
    if ((columns != 0 && rows > std::numeric_limits<IndexType>::max() / columns)
        || Arguments.lhs.size() != rows * columns || !IsValidTolerance(Arguments.options.dual_tolerance)
        || !IsValidTolerance(Arguments.options.rank_tolerance) || Overlaps(Arguments.lhs, Arguments.rhs)
        || Overlaps(Arguments.lhs, Arguments.solution) || Overlaps(Arguments.rhs, Arguments.solution)) {
        return { .status = Status::InvalidInput };
    }
    const auto is_finite = [](double Value) { return std::isfinite(Value); };
    if (!std::ranges::all_of(Arguments.lhs, is_finite) || !std::ranges::all_of(Arguments.rhs, is_finite)) {
        return { .status = Status::InvalidInput };
    }
    auto& r_workspace = Arguments.r_workspace;
    r_workspace.mDual.resize(columns);
    r_workspace.mColumnNorms.resize(columns);
    r_workspace.mTrial.resize(rows);
    r_workspace.mReflector.resize(rows);
    r_workspace.mIndices.resize(columns);
    return Solver(
               Arguments,
               r_workspace.mDual,
               r_workspace.mColumnNorms,
               r_workspace.mTrial,
               r_workspace.mReflector,
               r_workspace.mIndices
    )
        .Run();
}

}  // namespace queso::nnls
