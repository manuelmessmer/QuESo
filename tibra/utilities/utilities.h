// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef UTILITIES_H
#define UTILITIES_H

namespace tibra {

namespace utilities {

    /// @brief Swap function for two ptrs.
    /// @tparam T Type.
    /// @param a Pointer 1.
    /// @param b Pointer 2.
    template <typename T>
    static inline void swap(T& a, T& b) {
        T tmp = a;
        a = b;
        b = tmp;
    }

    /// @brief Swap function for two unique ptrs.
    /// @tparam T Type
    /// @param a Pointer 1.
    /// @param b Poitner 2.
    template <typename T>
    static inline void swap_unique(T& a, T& b) {
        T tmp = std::move(a);
        a = std::move(b);
        b = std::move(tmp);
    }

    ///@brief Simple Power functions
    ///@param x Value
    ///@param p Order
    ///@details For gcc (without --ffast-math compiler flag) this is faster than std::pow().
    static double inline power( double x, std::size_t p){
        double result = 1.0;
        while( p > 0UL ) {
            result = result * x;
            p -= 1;
        }
        return result;
    }

} // End namespace utiliites

} // End namespace tibra

#endif