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
#include <memory>
#include <numeric>
#include <vector>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/includes/define.hpp"
#include "queso/includes/dereference_iterator.hpp"

namespace queso::testing {

    TEST_SUITE("DereferenceIteratorTestSuite")
    {

        TEST_CASE("DereferenceIteratorTestSuite - Test Dereference Iterator Member Operations")
        {


            {  // Iterate
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(1));
                vec.emplace_back(MakeUnique<int>(2));
                vec.emplace_back(MakeUnique<int>(3));

                auto begin = dereference_iterator(vec.begin());
                auto end = dereference_iterator(vec.end());

                std::vector<int> values;
                for (auto it = begin; it != end; ++it) { values.push_back(*it); }

                REQUIRE((values.size()) == (3));
                CHECK((values[0]) == (1));
                CHECK((values[1]) == (2));
                CHECK((values[2]) == (3));
            }
            {  // operator= and operator*
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(1));
                vec.emplace_back(MakeUnique<int>(2));

                auto it1 = dereference_iterator(vec.begin());
                auto it2 = dereference_iterator(vec.end());

                it2 = it1;
                CHECK(it2 == it1);
                CHECK((*it2) == (1));
            }
            {  // operator[]
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(11));
                vec.emplace_back(MakeUnique<int>(22));
                vec.emplace_back(MakeUnique<int>(33));

                auto begin = dereference_iterator(vec.begin());

                CHECK((begin[0]) == (11));
                CHECK((begin[1]) == (22));
                CHECK((begin[2]) == (33));
            }
            {  // operator->
                struct Data
                {
                    int value;
                };
                std::vector<std::unique_ptr<Data>> vec;
                vec.emplace_back(MakeUnique<Data>(Data{ 42 }));
                vec.emplace_back(MakeUnique<Data>(Data{ 43 }));

                auto it = dereference_iterator(vec.begin());

                CHECK((it->value) == (42));
                ++it;
                CHECK((it->value) == (43));
            }
            {  // operator++ and operator--
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(100));
                vec.emplace_back(MakeUnique<int>(200));
                vec.emplace_back(MakeUnique<int>(300));

                auto it = dereference_iterator(vec.begin());
                auto end = dereference_iterator(vec.end());

                CHECK((*it) == (100));

                // Test pre-increment (++it)
                auto& pre_inc = ++it;
                CHECK((*pre_inc) == (200));
                CHECK(&pre_inc == &it);  // pre-increment returns reference to self

                // Test post-increment (it++)
                auto post_inc = it++;
                CHECK((*post_inc) == (200));  // post-increment returns value before increment
                CHECK((*it) == (300));
                CHECK(&post_inc != &it);  // post-increment returns a copy

                // Test pre-decrement (--it)
                auto& pre_dec = --it;
                CHECK((*pre_dec) == (200));
                CHECK(&pre_dec == &it);

                // Test post-decrement (it--)
                auto post_dec = it--;
                CHECK((*post_dec) == (200));
                CHECK((*it) == (100));
                CHECK(&post_dec != &it);

                CHECK(it != end);
            }
            {  // operator+= and operator-=
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(10));
                vec.emplace_back(MakeUnique<int>(20));
                vec.emplace_back(MakeUnique<int>(30));
                vec.emplace_back(MakeUnique<int>(40));

                auto it = dereference_iterator(vec.begin());
                auto end = dereference_iterator(vec.end());

                it += 2;
                CHECK((*it) == (30));
                it -= 1;
                CHECK((*it) == (20));
                it += 2;
                CHECK((*it) == (40));
                it -= 3;
                CHECK((*it) == (10));
                CHECK(it != end);
            }
        }

        TEST_CASE("DereferenceIteratorTestSuite - Test Dereference Iterator Non Member Operations")
        {


            {  // operator+ and operator-
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(5));
                vec.emplace_back(MakeUnique<int>(15));
                vec.emplace_back(MakeUnique<int>(25));

                auto it = dereference_iterator(vec.begin());

                auto it2 = it + 2;
                CHECK((*it2) == (25));

                auto it3 = it2 - 1;
                CHECK((*it3) == (15));

                CHECK((1 + it3) == it2);
                CHECK((it2 - 2) == it);
            }

            {  // operator- (with two iterators)
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(1));
                vec.emplace_back(MakeUnique<int>(2));
                vec.emplace_back(MakeUnique<int>(3));
                vec.emplace_back(MakeUnique<int>(4));

                auto begin = dereference_iterator(vec.begin());
                auto mid = begin + 2;
                auto end = dereference_iterator(vec.end());

                CHECK((end - begin) == (4));
                CHECK((mid - begin) == (2));
                CHECK((end - mid) == (2));
            }

            {  // operator<, operator<=, operator>, operator>=, operator==, operator!=
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(1));
                vec.emplace_back(MakeUnique<int>(2));
                vec.emplace_back(MakeUnique<int>(3));

                auto begin = dereference_iterator(vec.begin());
                auto mid = begin + 1;
                auto end = dereference_iterator(vec.end());

                // operator<
                CHECK(begin < mid);
                CHECK(begin < end);
                CHECK(!(mid < begin));
                CHECK(!(end < end));

                // operator<=
                CHECK(begin <= begin);
                CHECK(begin <= mid);
                CHECK(mid <= end);
                CHECK(end <= end);
                CHECK(!(mid <= begin));

                // operator>
                CHECK(mid > begin);
                CHECK(end > begin);
                CHECK(!(begin > mid));
                CHECK(!(begin > begin));

                // operator>=
                CHECK(mid >= begin);
                CHECK(end >= mid);
                CHECK(end >= end);
                CHECK(begin >= begin);
                CHECK(!(begin >= mid));

                // operator==
                CHECK(begin == begin);
                CHECK(!(begin == mid));
                CHECK(end == end);

                // operator!=
                CHECK(begin != mid);
                CHECK(!(begin != begin));
                CHECK(mid != end);
                CHECK(!(end != end));
            }
        }

        TEST_CASE("DereferenceIteratorTestSuite - Dereference Iterator With STLAlgorithms")
        {

            {  // std::find
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(10));
                vec.emplace_back(MakeUnique<int>(20));
                vec.emplace_back(MakeUnique<int>(30));

                auto begin = dereference_iterator(vec.begin());
                auto end = dereference_iterator(vec.end());

                auto it = std::find(begin, end, 20);
                CHECK(it != end);
                CHECK((*it) == (20));
            }

            {  // std::sort
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(30));
                vec.emplace_back(MakeUnique<int>(10));
                vec.emplace_back(MakeUnique<int>(20));

                auto begin = dereference_iterator(vec.begin());
                auto end = dereference_iterator(vec.end());

                // Sort in ascending order
                std::sort(begin, end);

                std::vector<int> sorted_values;
                for (auto it = begin; it != end; ++it) { sorted_values.push_back(*it); }

                REQUIRE((sorted_values.size()) == (3));
                CHECK((sorted_values[0]) == (10));
                CHECK((sorted_values[1]) == (20));
                CHECK((sorted_values[2]) == (30));
            }

            {  // std::copy
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(1));
                vec.emplace_back(MakeUnique<int>(2));
                vec.emplace_back(MakeUnique<int>(3));

                auto begin = dereference_iterator(vec.begin());
                auto end = dereference_iterator(vec.end());

                std::vector<int> out(3);
                std::copy(begin, end, out.begin());

                CHECK((out[0]) == (1));
                CHECK((out[1]) == (2));
                CHECK((out[2]) == (3));
            }

            {  // std::accumulate
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(4));
                vec.emplace_back(MakeUnique<int>(5));
                vec.emplace_back(MakeUnique<int>(6));

                auto begin = dereference_iterator(vec.begin());
                auto end = dereference_iterator(vec.end());

                int sum = std::accumulate(begin, end, 0);
                CHECK((sum) == (15));
            }

            {  // std::distance and std::advance
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(10));
                vec.emplace_back(MakeUnique<int>(20));
                vec.emplace_back(MakeUnique<int>(30));
                vec.emplace_back(MakeUnique<int>(40));

                auto begin = dereference_iterator(vec.begin());
                auto end = dereference_iterator(vec.end());

                auto dist = std::distance(begin, end);
                CHECK((dist) == (4));

                auto it = begin;
                std::advance(it, 2);
                CHECK((*it) == (30));
            }
        }

        TEST_CASE("DereferenceIteratorTestSuite - Test Dereference Iterator Const Correctness")
        {

            {
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(5));
                vec.emplace_back(MakeUnique<int>(6));

                auto begin = dereference_iterator(vec.begin());
                auto end = dereference_iterator(vec.end());
                CHECK(begin != end);

                auto& v1 = *begin;
                static_assert(
                    !std::is_const<std::remove_reference_t<decltype(v1)>>::value, "v must be a non-const reference"
                );
                auto& v2 = begin[0];
                static_assert(
                    !std::is_const<std::remove_reference_t<decltype(v2)>>::value, "v must be a non-const reference"
                );

                auto c_begin = dereference_iterator(vec.cbegin());
                auto c_end = dereference_iterator(vec.cend());

                auto& c_v1 = *c_begin;
                static_assert(
                    std::is_const<std::remove_reference_t<decltype(c_v1)>>::value, "v must be a const reference"
                );
                auto& c_v2 = c_begin[0];
                static_assert(
                    std::is_const<std::remove_reference_t<decltype(c_v2)>>::value, "v must be a const reference"
                );

                int sum = 0;
                for (auto it = c_begin; it != c_end; ++it) { sum += *it; }
                CHECK((sum) == (11));
            }
            {  // operator->
                struct Data
                {
                    int value;
                };
                std::vector<std::unique_ptr<Data>> vec;
                vec.emplace_back(MakeUnique<Data>(Data{ 42 }));
                vec.emplace_back(MakeUnique<Data>(Data{ 43 }));

                auto it = dereference_iterator(vec.begin());
                auto& v = it->value;
                static_assert(
                    !std::is_const<std::remove_reference_t<decltype(v)>>::value, "v must be a non-const reference"
                );

                auto c_it = dereference_iterator(vec.cbegin());
                auto& c_v = c_it->value;
                static_assert(
                    std::is_const<std::remove_reference_t<decltype(c_v)>>::value, "v must be a const reference"
                );
            }
        }

        TEST_CASE("DereferenceIteratorTestSuite - Test Dereference Range")
        {


            {  // Test DereferenceRange basic usage
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(100));
                vec.emplace_back(MakeUnique<int>(200));
                vec.emplace_back(MakeUnique<int>(300));

                auto range = dereference_range(vec.begin(), vec.end());
                std::vector<int> values;
                for (auto v : range) { values.push_back(v); }
                REQUIRE((values.size()) == (3));
                CHECK((values[0]) == (100));
                CHECK((values[1]) == (200));
                CHECK((values[2]) == (300));
            }

            {  // Test DereferenceRange with const_iterator
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(1));
                vec.emplace_back(MakeUnique<int>(2));

                const auto& cvec = vec;
                auto range = dereference_range(cvec.cbegin(), cvec.cend());
                int sum = 0;
                for (auto& v : range) {
                    static_assert(
                        std::is_const<std::remove_reference_t<decltype(v)>>::value, "v must be a const reference"
                    );
                    sum += v;
                }
                CHECK((sum) == (3));
            }

            {  // Test DereferenceRange with empty range
                std::vector<std::unique_ptr<int>> vec;
                auto range = dereference_range(vec.begin(), vec.end());
                int count = 0;
                for (auto v : range) {
                    (void)v;
                    ++count;
                }
                CHECK((count) == (0));
            }

            {  // Test DereferenceRange with STL algorithms (std::accumulate)
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(10));
                vec.emplace_back(MakeUnique<int>(20));
                vec.emplace_back(MakeUnique<int>(30));

                auto range = dereference_range(vec.begin(), vec.end());
                int sum = std::accumulate(range.begin(), range.end(), 0);
                CHECK((sum) == (60));
            }

            {  // Test DereferenceRange with std::find
                std::vector<std::unique_ptr<int>> vec;
                vec.emplace_back(MakeUnique<int>(7));
                vec.emplace_back(MakeUnique<int>(8));
                vec.emplace_back(MakeUnique<int>(9));

                auto range = dereference_range(vec.begin(), vec.end());
                auto it = std::find(range.begin(), range.end(), 8);
                CHECK(it != range.end());
                CHECK((*it) == (8));
            }
        }

    }  // TEST_SUITE

} // namespace queso::testing
