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
#include <exception>
#include <mutex>
#include <utility>

namespace queso::detail {

/// @brief Captures a worker failure for rethrow after an OpenMP region.
/// @details Workers finish their remaining iterations. Results must not be consumed before Rethrow().
class ParallelException
{
public:
    ///@name Operations
    ///@{

    /// @brief Executes one iteration without allowing its exception to escape the worker.
    /// @tparam TOperation Nullary callable type.
    /// @param rOperation Iteration body, invoked once with its original value category; partially modified results
    ///                   are discarded by the caller on failure.
    template<typename TOperation>
    void Execute(TOperation&& rOperation)
    {
        try {
            std::forward<TOperation>(rOperation)();
        } catch (...) {
            const std::lock_guard lock(mMutex);
            if (!mException) { mException = std::current_exception(); }
        }
    }

    /// @brief Rethrows a captured failure. Call only after all workers have joined.
    void Rethrow() const
    {
        if (mException) { std::rethrow_exception(mException); }
    }

    ///@}
private:
    ///@name Private Members
    ///@{
    std::mutex mMutex;
    std::exception_ptr mException;
    ///@}
};

}  // namespace queso::detail
