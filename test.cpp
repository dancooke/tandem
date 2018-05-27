/*  tandem.hpp
 
 Copyright (C) 2017-2018 University of Oxford.
 
 Author: Daniel Cooke <dcooke@well.ox.ac.uk>
 
 Use of this source code is governed by the MIT license that can be found in the LICENSE file. */


#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>

#include "tandem.hpp"

template <typename Container1, typename Container2>
bool are_equal(const Container1& lhs, const Container2& rhs)
{
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::cend(rhs));
}

auto get_substrings(const std::string& str, const std::vector<tandem::LZBlock>& lz_blocks)
{
    std::vector<std::string> result {};
    result.reserve(lz_blocks.size());
    std::transform(std::cbegin(lz_blocks), std::cend(lz_blocks), std::back_inserter(result),
                   [&str] (const auto& block) { return str.substr(block.pos, block.length); });
    return result;
}

auto lz_decomposition(std::string str)
{
    str.push_back('$'); // sentinel
    auto result = get_substrings(str, tandem::lempel_ziv_factorisation(str));
    assert(!result.empty() && result.back() == "$");
    result.pop_back();
    return result;
}

int test_lz_decomposition()
{
    int num_failed_tests {0};
    
    const std::string seq1 {"abcdbabdcadbcbdbabcabcb"};
    const auto lz_factorization1 = lz_decomposition(seq1);
    const std::vector<std::string> true_lz_factorization1 {"a", "b", "c", "d", "b", "ab", "d", "c", "a", "db", "c", "bd", "bab", "ca", "bcb"};
    if (!are_equal(lz_factorization1, true_lz_factorization1)) {
        ++num_failed_tests;
    }
    
    const std::string seq2 {"aaaaa"};
    const auto lz_factorization2 = lz_decomposition(seq2);
    const std::vector<std::string> true_lz_factorization2 {"a", "aaaa"};
    if (!are_equal(lz_factorization2, true_lz_factorization2)) {
        ++num_failed_tests;
    }
    
    const std::string seq3 {"abcbcbbcbaba"};
    const auto lz_factorization3 = lz_decomposition(seq3);
    const std::vector<std::string> true_lz_factorization3 {"a", "b", "c", "bcb", "bcb", "ab", "a"};
    if (!are_equal(lz_factorization3, true_lz_factorization3)) {
        ++num_failed_tests;
    }
    
    const std::string seq4 {"AAATGAATT"};
    const auto lz_factorization4 = lz_decomposition(seq4);
    const std::vector<std::string> true_lz_factorization4 {"A", "AA", "T", "G", "AAT", "T"};
    if (!are_equal(lz_factorization4, true_lz_factorization4)) {
        ++num_failed_tests;
    }
    
    return num_failed_tests;
}

auto lz_decomposition(std::string str, std::vector<std::uint32_t>& prev_occs)
{
    str.push_back('$'); // sentinel
    std::vector<tandem::LZBlock> lz_blocks;
    std::tie(lz_blocks, prev_occs) = tandem::lempel_ziv_factorisation_with_prev_block_occurences(str);
    auto result = get_substrings(str, lz_blocks);
    assert(!result.empty() && result.back() == "$");
    result.pop_back();
    return result;
}

int test_lz_decomposition_with_prev_blocks()
{
    int num_failed_tests {0};
    
    std::vector<std::uint32_t> prev_occs;
    
    const std::string seq1 {"abcdbabdcadbcbdbabcabcb"};
    const auto lz_factorization1 = lz_decomposition(seq1, prev_occs);
    const std::vector<std::string> true_lz_factorization1 {"a", "b", "c", "d", "b", "ab", "d", "c", "a", "db", "c", "bd", "bab", "ca", "bcb"};
    if (!are_equal(lz_factorization1, true_lz_factorization1)) {
        ++num_failed_tests;
    }
    
    const std::string seq2 {"aaaaa"};
    const auto lz_factorization2 = lz_decomposition(seq2, prev_occs);
    const std::vector<std::string> true_lz_factorization2 {"a", "aaaa"};
    if (!are_equal(lz_factorization2, true_lz_factorization2)) {
        ++num_failed_tests;
    }
    
    const std::string seq3 {"abcbcbbcbaba"};
    const auto lz_factorization3 = lz_decomposition(seq3, prev_occs);
    const std::vector<std::string> true_lz_factorization3 {"a", "b", "c", "bcb", "bcb", "ab", "a"};
    if (!are_equal(lz_factorization3, true_lz_factorization3)) {
        ++num_failed_tests;
    }
    
    const std::string seq4 {"AAATGAAT"};
    const auto lz_factorization4 = lz_decomposition(seq4, prev_occs);
    const std::vector<std::string> true_lz_factorization4 {"A", "AA", "T", "G", "AAT"};
    if (!are_equal(lz_factorization4, true_lz_factorization4)) {
        ++num_failed_tests;
    }
    
    return num_failed_tests;
}

auto get_mains_leftmost_maximal_periodiciites(std::string str)
{
    str.push_back('$');
    const auto lz_blocks = tandem::lempel_ziv_factorisation(str);
    auto lmrs = tandem::detail::find_leftmost_maximal_repetitions(str, lz_blocks);
    assert(!lmrs.empty());
    std::vector<std::string> result {};
    result.reserve(lmrs.size());
    for (const auto& lmr : lmrs) {
        result.push_back(str.substr(lmr.pos, lmr.length));
    }
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

int test_mains_algorithm()
{
    int num_failed_tests {0};
    
    const std::string seq1 {"abcbcbbcbaba"};
    const auto lmrs1 = get_mains_leftmost_maximal_periodiciites(seq1);
    const std::vector<std::string> true_lmrs1 {"baba", "bb", "bcbbcb", "bcbcb"};
    if (!are_equal(lmrs1, true_lmrs1)) {
        ++num_failed_tests;
    }
    
    const std::string seq2 {"AAATGAAT"};
    const auto lmrs2 = get_mains_leftmost_maximal_periodiciites(seq2);
    const std::vector<std::string> true_lmrs2 {"AA", "AAA"};
    if (!are_equal(lmrs2, true_lmrs2)) {
        ++num_failed_tests;
    }
    
    return num_failed_tests;
}

auto get_repeat_strings(std::string str)
{
    str.push_back('$');
    auto repeats = tandem::extract_exact_tandem_repeats(str);
    std::vector<std::string> result {};
    result.reserve(repeats.size());
    for (const auto& repeat : repeats) {
        result.push_back(str.substr(repeat.pos, repeat.length));
    }
    return result;
}

int test_tandem_repeat_finder()
{
    int num_failed_tests {0};
    
    const std::string seq1 {"AAATGAAT$"};
    const auto repeats1 = get_repeat_strings(seq1);
    const std::vector<std::string> true_repeats1 {"AAA", "AA"};
    if (are_equal(repeats1, true_repeats1)) {
        ++num_failed_tests;
    }
    
    std::ifstream file {"../test_dna.txt"};
    assert(file.good());
    std::string seq2;
    std::getline(file, seq2);
    seq2.push_back('$');
    const auto repeats2 = tandem::extract_exact_tandem_repeats(seq2);
    for (const auto& repeat : repeats2) {
        const auto substr = seq2.substr(repeat.pos, repeat.length);
        if (!std::equal(substr.cbegin(), substr.cbegin() + repeat.period, substr.cbegin() + repeat.period)) {
            ++num_failed_tests;
            break;
        }
    }
    
    return num_failed_tests;
}

 int main()
{
    int num_failed_tests {0};
    
    num_failed_tests += test_lz_decomposition();
    num_failed_tests += test_lz_decomposition_with_prev_blocks();
    num_failed_tests += test_mains_algorithm();
    num_failed_tests += test_tandem_repeat_finder();
    
    if (num_failed_tests > 0) {
        std::cout << num_failed_tests << " tests FAILED!" << std::endl;
    } else {
        std::cout << "All tests PASSED!" << std::endl;
    }
    
}
