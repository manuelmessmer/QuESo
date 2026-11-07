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

#pragma once

//// STL includes
#include <limits>
#include <span>
#include <vector>

//// Project includes
#include "queso/includes/core_definitions.hpp"

namespace queso::nnls {

/// @brief Controls the iteration budget and relative numerical tests. Zero selects automatic defaults.
struct Options
{
    /// Maximum number of active triangular solves. Automatic is 5 * number of columns, saturated on overflow.
    IndexType max_iterations = 0;
    /// Relative dual threshold: w_j / (||a_j||_2 * ||residual||_2). Automatic is 8 * epsilon * max(1, rows).
    double dual_tolerance = 0.0;
    /// Relative independence threshold: trailing column norm / original column norm. Automatic is epsilon * max(1,
    /// rows).
    double rank_tolerance = 0.0;
};

/// @brief Distinguishes optimal solutions from incomplete solves and invalid problems.
enum class Status {
    Converged,
    IterationLimit,
    InvalidInput,
    NumericalFailure,
};

/// @brief Outcome of an NNLS solve. Converged and IterationLimit guarantee a usable feasible solution.
struct Result
{
    Status status = Status::InvalidInput;
    /// Full residual of the returned iterate, including coherent numerical failures. Infinity means unavailable.
    double residual_norm = std::numeric_limits<double>::infinity();
    /// Number of attempted active triangular solves, including feasibility-restoration solves.
    IndexType iterations = 0;
};

struct Parameters;

/// @brief Caller-owned O(rows + columns) scratch storage, reusable across changing problem sizes.
/// @details Capacity is retained. Concurrent solves must use distinct workspaces.
class Workspace
{
private:
    friend Result Solve(Parameters Arguments);

    ///@name Private Members
    ///@{
    std::vector<double> mDual;
    std::vector<double> mColumnNorms;
    std::vector<double> mTrial;
    std::vector<double> mReflector;
    std::vector<IndexType> mIndices;
    ///@}
};

/// @brief Non-owning, disjoint problem buffers and a reusable workspace for an in-place cold-start solve.
struct Parameters
{
    /// Column-major matrix, overwritten with orthogonally transformed columns. Its size must equal rhs.size() *
    /// solution.size().
    std::span<double> lhs;
    /// Right-hand side, overwritten with its orthogonal transform. Its size defines the number of rows.
    std::span<double> rhs;
    /// Pre-sized output defining the number of columns. Initial values are ignored. Inactive entries are exactly zero.
    std::span<double> solution;
    Workspace& r_workspace;
    Options options{};
};

/// @brief Minimizes ||lhs * x - rhs||_2 subject to x >= 0 using Lawson-Hanson incremental QR active-set updates.
/// @details Overwrites lhs and rhs. Buffers must remain valid and disjoint. Invalid input leaves them unchanged.
/// Nonzero subnormal column or normalization norms return NumericalFailure. Trivial empty or zero-RHS solves are
/// exempt. Tolerances must be finite and in [0, 1). Empty equation sets return a zero solution. Empty variable sets
/// return
/// ||rhs||_2. Reuses workspace capacity. Allocation failures may throw. Numerical failure may leave transformed inputs.
/// @param Arguments Problem buffers, workspace and optional controls. Spans do not own their storage.
/// @return Status, actual residual norm when usable, and consumed iteration count. No logging is performed.
[[nodiscard]] Result Solve(Parameters Arguments);

}  // namespace queso::nnls
