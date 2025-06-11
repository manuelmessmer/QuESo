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

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
#include <numeric>

//// Project includes
#include "queso/includes/checks.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( DereferenceIteratorTestSuite )

BOOST_AUTO_TEST_CASE(TestDereferenceIteratorMemberOperations) {
    QuESo_INFO << "Testing :: Test DereferenceIterator :: Member Operations" << std::endl;

    { // Iterate
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(1));
        vec.emplace_back(MakeUnique<int>(2));
        vec.emplace_back(MakeUnique<int>(3));

        auto begin = dereference_iterator(vec.begin());
        auto end = dereference_iterator(vec.end());

        std::vector<int> values;
        for (auto it = begin; it != end; ++it) {
            values.push_back(*it);
        }

        QuESo_CHECK_EQUAL(values.size(), 3);
        QuESo_CHECK_EQUAL(values[0], 1);
        QuESo_CHECK_EQUAL(values[1], 2);
        QuESo_CHECK_EQUAL(values[2], 3);
    }
    { // operator= and operator*
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(1));
        vec.emplace_back(MakeUnique<int>(2));

        auto it1 = dereference_iterator(vec.begin());
        auto it2 = dereference_iterator(vec.end());

        it2 = it1;
        QuESo_CHECK(it2 == it1);
        QuESo_CHECK_EQUAL(*it2, 1);
    }
    { // operator[]
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(11));
        vec.emplace_back(MakeUnique<int>(22));
        vec.emplace_back(MakeUnique<int>(33));

        auto begin = dereference_iterator(vec.begin());

        QuESo_CHECK_EQUAL(begin[0], 11);
        QuESo_CHECK_EQUAL(begin[1], 22);
        QuESo_CHECK_EQUAL(begin[2], 33);
    }
    { // operator->
        struct Data { int value; };
        std::vector<std::unique_ptr<Data>> vec;
        vec.emplace_back(MakeUnique<Data>(Data{42}));
        vec.emplace_back(MakeUnique<Data>(Data{43}));

        auto it = dereference_iterator(vec.begin());

        QuESo_CHECK_EQUAL(it->value, 42);
        ++it;
        QuESo_CHECK_EQUAL(it->value, 43);
    }
    { // operator++ and operator--
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(100));
        vec.emplace_back(MakeUnique<int>(200));
        vec.emplace_back(MakeUnique<int>(300));

        auto it = dereference_iterator(vec.begin());
        auto end = dereference_iterator(vec.end());

        QuESo_CHECK_EQUAL(*it, 100);

        // Test pre-increment (++it)
        auto& pre_inc = ++it;
        QuESo_CHECK_EQUAL(*pre_inc, 200);
        QuESo_CHECK(&pre_inc == &it); // pre-increment returns reference to self

        // Test post-increment (it++)
        auto post_inc = it++;
        QuESo_CHECK_EQUAL(*post_inc, 200); // post-increment returns value before increment
        QuESo_CHECK_EQUAL(*it, 300);
        QuESo_CHECK(&post_inc != &it); // post-increment returns a copy

        // Test pre-decrement (--it)
        auto& pre_dec = --it;
        QuESo_CHECK_EQUAL(*pre_dec, 200);
        QuESo_CHECK(&pre_dec == &it);

        // Test post-decrement (it--)
        auto post_dec = it--;
        QuESo_CHECK_EQUAL(*post_dec, 200);
        QuESo_CHECK_EQUAL(*it, 100);
        QuESo_CHECK(&post_dec != &it);

        QuESo_CHECK(it != end);
    }
    { // operator+= and operator-=
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(10));
        vec.emplace_back(MakeUnique<int>(20));
        vec.emplace_back(MakeUnique<int>(30));
        vec.emplace_back(MakeUnique<int>(40));

        auto it = dereference_iterator(vec.begin());
        auto end = dereference_iterator(vec.end());

        it += 2;
        QuESo_CHECK_EQUAL(*it, 30);
        it -= 1;
        QuESo_CHECK_EQUAL(*it, 20);
        it += 2;
        QuESo_CHECK_EQUAL(*it, 40);
        it -= 3;
        QuESo_CHECK_EQUAL(*it, 10);
        QuESo_CHECK(it != end);
    }
}

BOOST_AUTO_TEST_CASE(TestDereferenceIteratorNonMemberOperations) {
    QuESo_INFO << "Testing :: Test DereferenceIterator :: Non-Member Operations" << std::endl;

    { // operator+ and operator-
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(5));
        vec.emplace_back(MakeUnique<int>(15));
        vec.emplace_back(MakeUnique<int>(25));

        auto it = dereference_iterator(vec.begin());

        auto it2 = it + 2;
        QuESo_CHECK_EQUAL(*it2, 25);

        auto it3 = it2 - 1;
        QuESo_CHECK_EQUAL(*it3, 15);

        QuESo_CHECK((1 + it3) == it2);
        QuESo_CHECK((it2 - 2) == it);
    }

    { // operator- (with two iterators)
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(1));
        vec.emplace_back(MakeUnique<int>(2));
        vec.emplace_back(MakeUnique<int>(3));
        vec.emplace_back(MakeUnique<int>(4));

        auto begin = dereference_iterator(vec.begin());
        auto mid = begin + 2;
        auto end = dereference_iterator(vec.end());

        QuESo_CHECK_EQUAL(end - begin, 4);
        QuESo_CHECK_EQUAL(mid - begin, 2);
        QuESo_CHECK_EQUAL(end - mid, 2);
    }

    { // operator<, operator<=, operator>, operator>=, operator==, operator!=
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(1));
        vec.emplace_back(MakeUnique<int>(2));
        vec.emplace_back(MakeUnique<int>(3));

        auto begin = dereference_iterator(vec.begin());
        auto mid = begin + 1;
        auto end = dereference_iterator(vec.end());

        // operator<
        QuESo_CHECK(begin < mid);
        QuESo_CHECK(begin < end);
        QuESo_CHECK(!(mid < begin));
        QuESo_CHECK(!(end < end));

        // operator<=
        QuESo_CHECK(begin <= begin);
        QuESo_CHECK(begin <= mid);
        QuESo_CHECK(mid <= end);
        QuESo_CHECK(end <= end);
        QuESo_CHECK(!(mid <= begin));

        // operator>
        QuESo_CHECK(mid > begin);
        QuESo_CHECK(end > begin);
        QuESo_CHECK(!(begin > mid));
        QuESo_CHECK(!(begin > begin));

        // operator>=
        QuESo_CHECK(mid >= begin);
        QuESo_CHECK(end >= mid);
        QuESo_CHECK(end >= end);
        QuESo_CHECK(begin >= begin);
        QuESo_CHECK(!(begin >= mid));

        // operator==
        QuESo_CHECK(begin == begin);
        QuESo_CHECK(!(begin == mid));
        QuESo_CHECK(end == end);

        // operator!=
        QuESo_CHECK(begin != mid);
        QuESo_CHECK(!(begin != begin));
        QuESo_CHECK(mid != end);
        QuESo_CHECK(!(end != end));
    }

}

BOOST_AUTO_TEST_CASE(DereferenceIteratorWithSTLAlgorithms) {
    QuESo_INFO << "Testing :: Test DereferenceIterator :: STL Algorithms" << std::endl;
    { // std::find
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(10));
        vec.emplace_back(MakeUnique<int>(20));
        vec.emplace_back(MakeUnique<int>(30));

        auto begin = dereference_iterator(vec.begin());
        auto end = dereference_iterator(vec.end());

        auto it = std::find(begin, end, 20);
        QuESo_CHECK(it != end);
        QuESo_CHECK_EQUAL(*it, 20);
    }

    { // std::sort
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(30));
        vec.emplace_back(MakeUnique<int>(10));
        vec.emplace_back(MakeUnique<int>(20));

        auto begin = dereference_iterator(vec.begin());
        auto end = dereference_iterator(vec.end());

        // Sort in ascending order
        std::sort(begin, end);

        std::vector<int> sorted_values;
        for (auto it = begin; it != end; ++it) {
            sorted_values.push_back(*it);
        }

        QuESo_CHECK_EQUAL(sorted_values.size(), 3);
        QuESo_CHECK_EQUAL(sorted_values[0], 10);
        QuESo_CHECK_EQUAL(sorted_values[1], 20);
        QuESo_CHECK_EQUAL(sorted_values[2], 30);
    }

    { // std::copy
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(1));
        vec.emplace_back(MakeUnique<int>(2));
        vec.emplace_back(MakeUnique<int>(3));

        auto begin = dereference_iterator(vec.begin());
        auto end = dereference_iterator(vec.end());

        std::vector<int> out(3);
        std::copy(begin, end, out.begin());

        QuESo_CHECK_EQUAL(out[0], 1);
        QuESo_CHECK_EQUAL(out[1], 2);
        QuESo_CHECK_EQUAL(out[2], 3);
    }

    { // std::accumulate
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(4));
        vec.emplace_back(MakeUnique<int>(5));
        vec.emplace_back(MakeUnique<int>(6));

        auto begin = dereference_iterator(vec.begin());
        auto end = dereference_iterator(vec.end());

        int sum = std::accumulate(begin, end, 0);
        QuESo_CHECK_EQUAL(sum, 15);
    }

    { // std::distance and std::advance
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(10));
        vec.emplace_back(MakeUnique<int>(20));
        vec.emplace_back(MakeUnique<int>(30));
        vec.emplace_back(MakeUnique<int>(40));

        auto begin = dereference_iterator(vec.begin());
        auto end = dereference_iterator(vec.end());

        auto dist = std::distance(begin, end);
        QuESo_CHECK_EQUAL(dist, 4);

        auto it = begin;
        std::advance(it, 2);
        QuESo_CHECK_EQUAL(*it, 30);
    }
}

BOOST_AUTO_TEST_CASE(TestDereferenceIteratorConstCorrectness) {
    QuESo_INFO << "Testing :: Test DereferenceIterator :: Const Correctness" << std::endl;
    {
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(5));
        vec.emplace_back(MakeUnique<int>(6));

        auto begin = dereference_iterator(vec.begin());
        auto end = dereference_iterator(vec.end());
        QuESo_CHECK(begin != end);

        auto& v1 = *begin;
        static_assert(!std::is_const<std::remove_reference_t<decltype(v1)>>::value, "v must be a non-const reference");
        auto& v2 = begin[0];
        static_assert(!std::is_const<std::remove_reference_t<decltype(v2)>>::value, "v must be a non-const reference");

        auto c_begin = dereference_iterator(vec.cbegin());
        auto c_end = dereference_iterator(vec.cend());

        auto& c_v1 = *c_begin;
        static_assert(std::is_const<std::remove_reference_t<decltype(c_v1)>>::value, "v must be a const reference");
        auto& c_v2 = c_begin[0];
        static_assert(std::is_const<std::remove_reference_t<decltype(c_v2)>>::value, "v must be a const reference");

        int sum = 0;
        for (auto it = c_begin; it != c_end; ++it) {
            sum += *it;
        }
        QuESo_CHECK_EQUAL(sum, 11);
    }
    { // operator->
        struct Data { int value; };
        std::vector<std::unique_ptr<Data>> vec;
        vec.emplace_back(MakeUnique<Data>(Data{42}));
        vec.emplace_back(MakeUnique<Data>(Data{43}));

        auto it = dereference_iterator(vec.begin());
        auto& v = it->value;
        static_assert(!std::is_const<std::remove_reference_t<decltype(v)>>::value, "v must be a non-const reference");

        auto c_it = dereference_iterator(vec.cbegin());
        auto& c_v = c_it->value;
        static_assert(std::is_const<std::remove_reference_t<decltype(c_v)>>::value, "v must be a const reference");
    }
}

BOOST_AUTO_TEST_CASE(TestDereferenceRange) {
    QuESo_INFO << "Testing :: Test DereferenceRange :: Basic Usage and STL algorithms" << std::endl;

    { // Test DereferenceRange basic usage
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(100));
        vec.emplace_back(MakeUnique<int>(200));
        vec.emplace_back(MakeUnique<int>(300));

        auto range = dereference_range(vec.begin(), vec.end());
        std::vector<int> values;
        for (auto v : range) {
            values.push_back(v);
        }
        QuESo_CHECK_EQUAL(values.size(), 3);
        QuESo_CHECK_EQUAL(values[0], 100);
        QuESo_CHECK_EQUAL(values[1], 200);
        QuESo_CHECK_EQUAL(values[2], 300);
    }

    { // Test DereferenceRange with const_iterator
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(1));
        vec.emplace_back(MakeUnique<int>(2));

        const auto& cvec = vec;
        auto range = dereference_range(cvec.cbegin(), cvec.cend());
        int sum = 0;
        for (auto& v : range) {
            static_assert(std::is_const<std::remove_reference_t<decltype(v)>>::value, "v must be a const reference");
            sum += v;
        }
        QuESo_CHECK_EQUAL(sum, 3);
    }

    { // Test DereferenceRange with empty range
        std::vector<std::unique_ptr<int>> vec;
        auto range = dereference_range(vec.begin(), vec.end());
        int count = 0;
        for (auto v : range) {
            (void)v;
            ++count;
        }
        QuESo_CHECK_EQUAL(count, 0);
    }

    { // Test DereferenceRange with STL algorithms (std::accumulate)
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(10));
        vec.emplace_back(MakeUnique<int>(20));
        vec.emplace_back(MakeUnique<int>(30));

        auto range = dereference_range(vec.begin(), vec.end());
        int sum = std::accumulate(range.begin(), range.end(), 0);
        QuESo_CHECK_EQUAL(sum, 60);
    }

    { // Test DereferenceRange with std::find
        std::vector<std::unique_ptr<int>> vec;
        vec.emplace_back(MakeUnique<int>(7));
        vec.emplace_back(MakeUnique<int>(8));
        vec.emplace_back(MakeUnique<int>(9));

        auto range = dereference_range(vec.begin(), vec.end());
        auto it = std::find(range.begin(), range.end(), 8);
        QuESo_CHECK(it != range.end());
        QuESo_CHECK_EQUAL(*it, 8);
    }
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso