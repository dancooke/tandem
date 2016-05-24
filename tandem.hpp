/*  tandem.hpp -- Finding tadem repeats in sequences
 
 Copyright (C) 2015 University of Oxford.
 
 Author: Daniel Cooke <dcooke@well.ox.ac.uk>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.  */

#ifndef __tandem__tandem__
#define __tandem__tandem__

#include <vector>
#include <deque>
#include <map>
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <utility>
#include <numeric>

#include <cassert>

#include "libdivsufsort/divsufsort.h"

/**
 This is a small library for finding exact tandem repeats in sequences. It also provides a
 simple API for making common string processing data structures, that are used in the repeat
 finding algorithm. The main method is:
 
 template <typename T>
 std::vector<StringRun> find_maximal_repetitions(const T& str, uint32_t min_period, uint32_t max_period)
 
 where T can be any contiguous sequence container with value_type char (e.g. std::string, std::vector<char>)
 that has a size() member function (so no const char*).
 
 The library only works with uint32_t at the moment. This is mostly because divsufsort does not use templates,
 so I'd either have to use uint32_t or uint64_t (or size_t).
 */

namespace Tandem
{
struct StringRun
{
    StringRun() = default;
    explicit StringRun(uint32_t pos, uint32_t length, uint32_t period)
    : pos {pos}, length {length}, period {period} {}
    
    uint32_t pos, length, period;
};

inline bool operator==(const StringRun& lhs, const StringRun& rhs)
{
    return lhs.pos == rhs.pos && lhs.length == rhs.length;
}
inline bool operator!=(const StringRun& lhs, const StringRun& rhs)
{
    return !operator==(lhs, rhs);
}

// Wrapper for divsufsort
template <typename T>
std::vector<uint32_t>& make_suffix_array(const T& str, std::vector<uint32_t>& result)
{
    std::vector<saidx_t> sa(str.size());
    
    divsufsort(str.data(), sa.data(), static_cast<int>(str.size()));
    
    result.resize(str.size());
    
    std::copy(std::cbegin(sa), std::cend(sa), std::begin(result));
    
    return result;
}

template <typename T>
std::vector<uint32_t> make_suffix_array(const T& str)
{
    std::vector<uint32_t> result(str.size());
    return make_suffix_array(str, result);
}

// rank array is inverse suffix array
std::vector<uint32_t>& make_rank_array(const std::vector<uint32_t>& suffix_array,
                                       std::vector<uint32_t>& result);
std::vector<uint32_t> make_rank_array(const std::vector<uint32_t>& suffix_array);

namespace detail
{
    // Llie et al (2010) show that this direct computation of LCE actually outperforms other O(n) and O(log n)
    // methods in practice.
    
    template <typename T>
    auto forward_lce(const T& str, uint32_t i, uint32_t j, uint32_t n)
    {
        using std::cbegin; using std::distance; using std::next; using std::mismatch;
        const auto it = cbegin(str);
        return static_cast<uint32_t>(distance(next(it, i), mismatch(next(it, i), next(it, n), next(it, j)).first));
    }
    
    template <typename T>
    auto forward_lce(const T& str, uint32_t i, uint32_t j)
    {
        return forward_lce(str, i, j, static_cast<uint32_t>(str.size()));
    }
    
    template <typename T>
    auto backward_lce(const T& str, uint32_t i, uint32_t j, uint32_t n)
    {
        using std::crend; using std::distance; using std::prev;  using std::mismatch;
        const auto it = crend(str);
        return static_cast<uint32_t>(distance(prev(it, i + 1), mismatch(prev(it, i + 1), prev(it, n), prev(it, j + 1)).first));
    }
    
    template <typename T>
    auto backward_lce(const T& str, uint32_t i, uint32_t j)
    {
        return backward_lce(str, i, j, uint32_t {});
    }
} // namespace detail

// LCP = Longest Common Prefix. O(n) implementation given in Kasai et al (2001).
template <typename T>
std::vector<uint32_t>
make_lcp_array(const T& str, const std::vector<uint32_t>& suffix_array,
               std::vector<uint32_t>& result)
{
    const auto rank = make_rank_array(suffix_array);
    
    const auto str_len = str.size();
    
    assert(rank.size() == str_len);
    
    result.resize(str_len);
    
    for (uint32_t i {0}, h {0}; i < str_len; ++i) {
        if (rank[i] > 0) {
            h += detail::forward_lce(str, i + h, suffix_array[rank[i] - 1] + h);
            result[rank[i]] = h;
            if (h > 0) --h;
        }
    }
    
    return result;
}

template <typename T>
std::vector<uint32_t>
make_lcp_array(const T& str, const std::vector<uint32_t>& suffix_array)
{
    std::vector<uint32_t> result(str.size());
    return make_lcp_array(str, suffix_array, result);
}

// LPF = Longest Previous Factor
std::vector<uint32_t>
make_lpf_array(std::vector<uint32_t>& sa, std::vector<uint32_t>& lcp);

std::pair<std::vector<uint32_t>, std::vector<uint32_t>>
make_lpf_and_prev_occ_arrays(std::vector<uint32_t>& sa, std::vector<uint32_t>& lcp);

namespace detail
{
    template <typename T>
    void make_sa_and_lcp(const T& str, std::vector<uint32_t>& sa, std::vector<uint32_t>& lcp)
    {
        const auto n = str.size() + 1; // + 1 to avoid reallocation downstream
        sa.reserve(n);
        make_suffix_array(str, sa);
        lcp.reserve(n);
        make_lcp_array(str, sa, lcp);
    }
} // namespace detail

template <typename T>
auto make_lpf_array(const T& str)
{
    std::vector<uint32_t> sa {}, lcp {};
    detail::make_sa_and_lcp(str, sa, lcp);
    return make_lpf_array(sa, lcp);
}

template <typename T>
auto make_lpf_and_prev_occ_arrays(const T& str)
{
    std::vector<uint32_t> sa {}, lcp {};
    detail::make_sa_and_lcp(str, sa, lcp);
    return make_lpf_and_prev_occ_arrays(sa, lcp);
}

struct LZBlock
{
    LZBlock() = default;
    explicit LZBlock(uint32_t pos, uint32_t length) : pos {pos}, length {length} {}
    uint32_t pos, length;
};

// Implementation of algorithm found in Crochemore et al. (2008)
template <typename T>
std::vector<LZBlock> lempel_ziv_factorisation(const T& str)
{
    if (str.empty()) return {};
    
    const auto lpf = make_lpf_array(str);
    
    const auto str_len = str.size();
    
    assert(lpf.size() == str_len);
    
    std::vector<LZBlock> result {};
    result.reserve(str_len); // max possible blocks
    
    assert(lpf.front() == 0);
    
    uint32_t end {1};  // start at 1 because the first element of lpf is always 0
    
    result.emplace_back(0, end);
    
    while (end < str_len) {
        const auto m = std::max(uint32_t {1}, lpf[end]);
        result.emplace_back(end, m);
        end += m;
    }
    
    result.shrink_to_fit();
    
    return result;
}

// Implementation of algorithm found in Crochemore et al. (2008)
template <typename T>
std::pair<std::vector<LZBlock>, std::vector<uint32_t>>
lempel_ziv_factorisation_with_prev_block_occurences(const T& str)
{
    if (str.empty()) return {{}, {}};
    
    std::vector<uint32_t> lpf, prev_occ;
    std::tie(lpf, prev_occ) = make_lpf_and_prev_occ_arrays(str);
    
    const auto str_len = str.size();
    
    assert(lpf.size() == str_len);
    
    std::vector<LZBlock> lz_blocks {};
    lz_blocks.reserve(str_len); // max possible blocks
    
    std::vector<uint32_t> prev_lz_block_occurrence {};
    prev_lz_block_occurrence.reserve(str_len);
    
    assert(lpf.front() == 0);
    
    uint32_t end {1}; // start at 1 because the first element of lpf is always 0
    
    lz_blocks.emplace_back(0, end);
    prev_lz_block_occurrence.emplace_back(-1);
    
    while (end < str_len) {
        const auto m = std::max(uint32_t {1}, lpf[end]);
        lz_blocks.emplace_back(end, m);
        prev_lz_block_occurrence.emplace_back(prev_occ[end]);
        end += m;
    }
    
    lz_blocks.shrink_to_fit();
    prev_lz_block_occurrence.shrink_to_fit();
    
    return std::make_pair(std::move(lz_blocks), std::move(prev_lz_block_occurrence));
}

namespace detail
{
    using RunQueue = std::deque<StringRun>;
    
    // Implements Mains algorithm found in Main (1989). Obscure notation as in paper.
    template <typename T>
    RunQueue
    find_leftmost_maximal_repetitions(const T& str, const std::vector<LZBlock>& lz_blocks,
                                      const uint32_t min_period = 1, const uint32_t max_period = -1)
    {
        RunQueue result {};
        
        for (uint32_t h {1}; h < lz_blocks.size(); ++h) {
            const auto u   = lz_blocks[h].pos;
            const auto n   = lz_blocks[h].length;
            const auto m   = std::min(u, 2 * lz_blocks[h - 1].length + n);
            const auto t   = u - m;
            const auto end = u + n;
            
            // rightmax periodicities
            for (uint32_t j {min_period}; j <= std::min(n, max_period); ++j) {
                const auto ls = detail::backward_lce(str, u - 1, u + j - 1, t);
                const auto lp = detail::forward_lce(str, u + j, u, end);
                
                if (ls > 0 && ls + lp >= j && j + lp < n) {
                    result.emplace_back(u - ls, j + lp + ls, j);
                }
            }
            
            // leftmax periodicities
            for (uint32_t j {min_period}; j < std::min(m, max_period); ++j) {
                const auto ls = detail::backward_lce(str, u - j - 1, u - 1, t);
                const auto lp = detail::forward_lce(str, u, u - j, end);
                
                if (ls + lp >= j) {
                    result.emplace_back(u - (ls + j), j + lp + ls, j);
                }
            }
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    using RunBuckets = std::vector<std::vector<StringRun>>;
    
    void init_end_buckets(std::size_t n, const RunQueue& lmrs, RunBuckets& result);
    
    template <typename T>
    RunBuckets extract_runs_end_buckets(const T& str, const std::vector<LZBlock>& lz_blocks,
                                        const uint32_t min_period, const uint32_t max_period)
    {
        const auto lmrs = find_leftmost_maximal_repetitions(str, lz_blocks, min_period, max_period);
        
        RunBuckets result;
        
        init_end_buckets(str.size(), lmrs, result);
        
        assert(result.size() == str.size());
        
        for (const auto& run : lmrs) {
            auto& bucket = result[run.pos + run.length - 1];
            
            const auto it = std::find(std::begin(bucket), std::end(bucket), run);
            
            if (it == std::end(bucket)) {
                bucket.push_back(run);
            } else if (it->period < run.period) {
                *it = run;
            }
        }
        
        return result;
    }
    
    void remove_empty_buckets(RunBuckets& buckets);
    
    void init_sorted_buckets(std::size_t n, const RunBuckets& end_buckets, RunBuckets& result);
    
    template <typename T>
    RunBuckets extract_runs_with_bucket_sort(const T& str, const std::vector<LZBlock>& lz_blocks,
                                             const uint32_t min_period, const uint32_t max_period)
    {
        auto end_buckets = extract_runs_end_buckets(str, lz_blocks, min_period, max_period);
        
        remove_empty_buckets(end_buckets);
        
        RunBuckets result;
        
        init_sorted_buckets(str.size(), end_buckets, result);
        
        for (auto& bucket : end_buckets) {
            for (const auto& run : bucket) {
                result[run.pos].push_back(run);
            }
            bucket.clear();
            bucket.shrink_to_fit();
        }
        
        return result;
    }
    
    template <typename T>
    RunBuckets extract_runs_with_inplace_sort(const T& str, const std::vector<LZBlock>& lz_blocks,
                                              const uint32_t min_period, const uint32_t max_period)
    {
        std::vector<StringRun> lmrs;
        
        {
            // switch from std::deque to std::vector
            auto tmp_lmrs = find_leftmost_maximal_repetitions(str, lz_blocks, min_period, max_period);
            lmrs.assign(std::begin(tmp_lmrs), std::end(tmp_lmrs));
        }
        
        std::sort(std::begin(lmrs), std::end(lmrs),
                  [] (const auto& lhs, const auto& rhs) {
                      if (lhs.pos != rhs.pos) return lhs.pos < rhs.pos;
                      if (lhs.length != rhs.length) return lhs.length < rhs.length;
                      return lhs.period > rhs.period;
                  });
        
        lmrs.erase(std::unique(std::begin(lmrs), std::end(lmrs)), std::end(lmrs));
        
        lmrs.shrink_to_fit();
        
        RunBuckets result(str.size(), RunBuckets::value_type {});
        
        for (auto it = std::cbegin(lmrs); it != std::cend(lmrs);) {
            const auto cur_pos = it->pos;
            
            const auto it2 = std::find_if_not(std::next(it), std::cend(lmrs),
                                              [cur_pos] (const auto& run) {
                                                  return run.pos == cur_pos;
                                              });
            
            result[cur_pos].assign(it, it2);
            
            it = it2;
        }
        
        return result;
    }
    
    template <typename T>
    RunBuckets extract_run_buckets(const T& str, const std::vector<LZBlock>& lz_blocks,
                                   const uint32_t min_period, const uint32_t max_period)
    {
        return extract_runs_with_bucket_sort(str, lz_blocks, min_period, max_period);
    }
    
    template <typename T>
    RunBuckets find_maximal_repetitions(const T& str,
                                        const std::vector<LZBlock>& lz_blocks,
                                        const std::vector<uint32_t>& prev_lz_block_occurrence,
                                        const uint32_t min_period, const uint32_t max_period)
    {
        auto sorted_buckets = extract_run_buckets(str, lz_blocks, min_period, max_period);
        
        for (uint32_t k {0}; k < lz_blocks.size(); ++k) {
            const auto& block = lz_blocks[k];
            
            const auto block_end = block.pos + block.length;
            const auto delta     = block.pos - ((prev_lz_block_occurrence[k] == -1) ? 0 : prev_lz_block_occurrence[k]);
            const auto v         = block_end - delta;
            
            for (auto j = block.pos; j < block_end; ++j) {
                const auto& target = sorted_buckets[j - delta];
                
                const auto last_target_itr = std::lower_bound(std::cbegin(target), std::cend(target), v,
                                                              [] (const auto& run, const auto val) {
                                                                  return run.pos + run.length < val;
                                                              });
                
                const auto num_targets = std::distance(std::cbegin(target), last_target_itr);
                
                if (num_targets > 0) {
                    std::vector<StringRun> shifted_targets(num_targets);
                    
                    std::transform(std::cbegin(target), last_target_itr, std::begin(shifted_targets),
                                   [delta] (const auto& run) {
                                       return StringRun {run.pos + delta, run.length, run.period};
                                   });
                    
                    sorted_buckets[j].insert(std::begin(sorted_buckets[j]),
                                             std::cbegin(shifted_targets),
                                             std::cend(shifted_targets));
                }
            }
        }
        
        return sorted_buckets;

    }
    
    // Implements the algorithm described in Kolpakov & Kucherov (1999)
    template <typename T>
    RunBuckets find_maximal_repetitions(const T& str, const uint32_t min_period, const uint32_t max_period)
    {
        std::vector<LZBlock> lz_blocks;
        std::vector<uint32_t> prev_lz_block_occurrence;
        
        std::tie(lz_blocks, prev_lz_block_occurrence) = lempel_ziv_factorisation_with_prev_block_occurences(str);
        
        return find_maximal_repetitions(str, lz_blocks, prev_lz_block_occurrence, min_period, max_period);
    }
    
    std::size_t count_runs(const RunBuckets& buckets);
    
    template <typename ForwardIt>
    std::vector<StringRun>
    find_homopolymers(const ForwardIt first, const ForwardIt last)
    {
        std::vector<StringRun> result {};
        
        if (first == last) return result;
        
        result.reserve(std::min(static_cast<std::size_t>(std::distance(first, last)), std::size_t {1024}));
        
        for (auto curr = first; curr != last; ) {
            const auto it = std::adjacent_find(curr, last);
            
            if (it == last) break;
            
            const auto base = *it;
            
            const auto it2 = std::find_if_not(std::next(it), last,
                                              [base] (const auto b) { return b == base; });
            
            result.emplace_back(static_cast<std::uint32_t>(std::distance(first, it)),
                                static_cast<std::uint32_t>(std::distance(it, it2)),
                                std::uint32_t {1});
            
            curr = it2;
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    template <std::size_t N, typename ForwardIt>
    std::vector<StringRun>
    find_exact_tandem_repeats(const ForwardIt first, const ForwardIt last)
    {
        std::vector<StringRun> result {};
        
        if (std::distance(first, last) < 2 * N) return result;
        
        auto it1 = std::adjacent_find(first, last, std::not_equal_to<void> {});
        
        if (it1 == last) return result;
        
        result.reserve(std::min(std::distance(first, last) / N, std::size_t {1024}));
        
        for (auto it2 = std::next(it1, N); it2 < last; ) {
            const auto p = std::mismatch(it2, last, it1);
            
            if (p.second >= it2) {
                result.emplace_back(static_cast<uint32_t>(std::distance(first, it1)),
                                    static_cast<uint32_t>(std::distance(it1, p.first)),
                                    N);
                it1 = p.second;
            } else {
                ++it1;
            }
            
            it1 = std::adjacent_find(it1, last, std::not_equal_to<void> {});
            
            if (it1 == last) break;
            
            it2 = std::next(it1, N);
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    template <typename T>
    auto find_homopolymers(const T& str)
    {
        return find_homopolymers(std::cbegin(str), std::cend(str));
    }
    
    template <typename T>
    auto find_exact_dinucleotide_tandem_repeats(const T& str)
    {
        return find_exact_tandem_repeats<2>(std::cbegin(str), std::cend(str));
    }
    
    template <typename T>
    auto find_exact_trinucleotide_tandem_repeats(const T& str)
    {
        return find_exact_tandem_repeats<3>(std::cbegin(str), std::cend(str));
    }
    
    template <typename Container1, typename Container2>
    void append(Container2&& src, Container1& dst)
    {
        const auto it = dst.insert(std::end(dst),
                                   std::make_move_iterator(std::begin(src)),
                                   std::make_move_iterator(std::end(src)));
        
        std::inplace_merge(std::begin(dst), it, std::end(dst),
                           [] (const StringRun& lhs, const StringRun& rhs) {
                               return lhs.pos < rhs.pos;
                           });
    }
} // namespace detail

/**
 This is the main function for finding exact repeats in a sequence.
 */
template <typename T>
std::vector<StringRun>
find_maximal_repetitions(const T& str, uint32_t min_period = 1, const uint32_t max_period = -1)
{
    if (min_period == 0) ++min_period;
    
    if (str.empty() || str.size() < min_period) return {};
    
    if (min_period > max_period) {
        throw std::domain_error {"find_maximal_repetitions: given unsatisfiable condition min_period > max_period"};
    }
    
    if (max_period <= 3) { // the naive algorithm is faster in these cases
        if (min_period == max_period) {
            switch(min_period) {
                case 1: return detail::find_homopolymers(str);
                case 2: return detail::find_exact_dinucleotide_tandem_repeats(str);
                case 3: return detail::find_exact_trinucleotide_tandem_repeats(str);
            }
        }
        
        if (min_period == 1) { // known max_period >= 2
            auto result = detail::find_homopolymers(str);
            
            detail::append(detail::find_exact_dinucleotide_tandem_repeats(str), result);
            
            if (max_period == 3) {
                detail::append(detail::find_exact_trinucleotide_tandem_repeats(str), result);
            }
            
            return result;
        } else { // min_period == 2 && max_period == 3
            auto result = detail::find_exact_dinucleotide_tandem_repeats(str);
            
            detail::append(detail::find_exact_trinucleotide_tandem_repeats(str), result);
            
            return result;
        }
    }
    
    auto sorted_buckets = detail::find_maximal_repetitions(str, min_period, max_period);
    
    std::vector<StringRun> result {};
    result.reserve(detail::count_runs(sorted_buckets));
    
    for (auto& bucket : sorted_buckets) {
        result.insert(std::end(result), std::cbegin(bucket), std::cend(bucket));
        bucket.clear();
        bucket.shrink_to_fit();
    }
    
    return result;
}

/**
 Replaces all contiguous sub-sequences of c with a single c, inplace, and returns a map of
 each c position in the new sequence, and how many c's have been removed up to the first non-c
 base past the position. This is just a helper that can speed up repetition finding if the sequence
 contains long runs of characters that are not of interest (e.g. unknwown base 'N' in DNA/RNA sequence).
 
 If this function is used, the output StringRun's will need to be rebased to get the correct positions
 in the original sequence. The function rebase does this transformation.
 
 Example:
 std::string str {"NNNACGTNNTGCNANNNN"};
 auto n_shift_map = colapse(str, 'N'); // str is now "NACGTNTGCNAN", n_shift_map contains (0, 2), (4, 3), (9, 6)
 */
template <typename SequenceType>
std::map<std::size_t, std::size_t> collapse(SequenceType& sequence, const char c)
{
    std::map<std::size_t, std::size_t> result {};
    
    const auto last = std::end(sequence);
    
    std::size_t position {0}, num_removed {0};
    
    for (auto first = std::begin(sequence); first != last;) {
        const auto it1 = std::adjacent_find(first, last,
                                            [c] (const char lhs, const char rhs) {
                                                return lhs == c && lhs == rhs;
                                            });
        
        if (it1 == last) break;
        
        const auto it2 = std::find_if_not(it1, last, [c] (const char b) { return b == c; });
        
        position    += std::distance(first, it1);
        num_removed += std::distance(it1, it2) - 1;
        
        result.emplace(position, num_removed);
        
        first = it2;
    }
    
    if (!result.empty()) {
        sequence.erase(std::unique(std::next(std::begin(sequence), std::cbegin(result)->first), last,
                                   [c] (const char lhs, const char rhs) {
                                       return lhs == c && lhs == rhs;
                                   }),
                       last);
    }
    
    return result;
}

void rebase(std::vector<Tandem::StringRun>& runs, const std::map<size_t, size_t>& shift_map);

} // namespace Tandem

#endif /* defined(__tandem__tandem__) */
