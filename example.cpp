#include <iostream>
#include <string>
#include <vector>
#include <cstdint>   // std::uint32_t;
#include <cstddef>   // std::size_t
#include <algorithm> // std::find_if_not, std::unique
#include <iterator>  // std::begin, std::end, std::distance, std::next
#include <map>

#include "tandem.h"

/**
    This is a simple example showing how to find all exact tandem repeats in a DNA sequence using
    the Tandem API.
 */

// Replaces all contiguous sub-sequences of N's with a single N, inplace, and returns a map of
// each N position in the new sequence, and how many N's have been removed up to the first non-N
// base past the position
template <typename SequenceType>
std::map<std::size_t, std::size_t> collapse_ns(SequenceType& sequence)
{
    std::map<std::size_t, std::size_t> result {};
    
    auto last = std::end(sequence);
    std::size_t position {}, num_removed {};
    
    for (auto first = std::begin(sequence); first != last;) {
        auto it1 = std::adjacent_find(first, last, [] (char lhs, char rhs) { return lhs == 'N' && lhs == rhs; });
        
        if (it1 == last) break;
        
        auto it2 = std::find_if_not(it1, last, [] (char base) { return base == 'N'; });
        
        position    += std::distance(first, it1);
        num_removed += std::distance(it1, it2) - 1;
        
        result.emplace(position, num_removed);
        
        first = it2;
    }
    
    if (!result.empty()) {
        sequence.erase(std::unique(std::next(std::begin(sequence), std::cbegin(result)->first), last,
                                   [] (char lhs, char rhs) { return lhs == 'N' && lhs == rhs; }), last);
    }
    
    return result;
}

void rebase(std::vector<Tandem::StringRun>& runs, const std::map<std::size_t, std::size_t>& shift_map)
{
    if (shift_map.empty()) return;
    
    auto shift_map_it = std::cbegin(shift_map);
    for (auto& run : runs) {
        while (std::next(shift_map_it)->first <= run.pos) ++shift_map_it;
        run.pos += static_cast<decltype(run.pos)>(shift_map_it->second);
    }
}

template <typename SequenceType>
bool is_dna(const SequenceType& sequence)
{
    return sequence.find_first_not_of("ACGTN") == SequenceType::npos;
}

template <typename SequenceType>
bool is_rna(const SequenceType& sequence)
{
    return sequence.find_first_not_of("ACGUN") == SequenceType::npos;
}

void find_big_repeats(std::string& sequence)
{
    std::cout << "Looking for repeats in " << sequence.size() << "bp" << std::endl;
    
    if (sequence.back() != 'N') {
        sequence.reserve(sequence.size() + 1);
        sequence.push_back('N'); // last character must not be in alphabet
    }
    
    auto n_shift_map = collapse_ns(sequence); // makes search a lot faster
    
    auto repeats = Tandem::find_maximal_repetitions(sequence);
    
    rebase(repeats, n_shift_map);
    
    std::cout << "Found " << repeats.size() << " repeats" << std::endl;
    
    if (!repeats.empty()) {
        auto it = std::max_element(std::cbegin(repeats), std::cend(repeats),
                                   [] (const auto& lhs, const auto& rhs) { return lhs.period < rhs.period; });
        
        std::cout << "biggest period is " << it->period << " at position " << it->pos << " with length " << it->length << std::endl;
        
        it = std::max_element(std::cbegin(repeats), std::cend(repeats),
                             [] (const auto& lhs, const auto& rhs) { return lhs.length < rhs.length; });
        
        std::cout << "longest run is " << it->length << " at position " << it->pos << " with period " << it->period << std::endl;
    }
}

int main()
{
    std::string sequence {};
    
    std::getline(std::cin, sequence);
    
    if (!is_dna(sequence) && !is_rna(sequence)) {
        std::cout << "this example is only for DNA or RNA sequences" << std::endl;
        exit(0);
    }
    
    find_big_repeats(sequence);
    
    return 0;
}
