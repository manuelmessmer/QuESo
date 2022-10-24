// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef UTILITIES_H
#define UTILITIES_H

namespace utilities {

    template <typename T>
    static inline void swap(T& a, T& b) {
        T tmp = a;
        a = b;
        b = tmp;
    }

    template <typename T>
    static inline void swap_unique(T& a, T& b) {
        T tmp = std::move(a);
        a = std::move(b);
        b = std::move(tmp);
    }

} // End namespace utiliites

#endif