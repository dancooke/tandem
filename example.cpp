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
    Simple example showing how to find all exact tandem repeats in a DNA sequence using
    the Tandem API.
 */

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
    if (sequence.empty()) return;
    
    std::cout << "Looking for repeats in " << sequence.size() << "bp" << std::endl;
    
    if (sequence.back() != 'N') {
        sequence.reserve(sequence.size() + 1);
        sequence.push_back('N'); // last character must not be in alphabet
    }
    
    auto n_shift_map = Tandem::collapse(sequence, 'N'); // makes search a lot faster
    
    auto repeats = Tandem::find_maximal_repetitions(sequence);
    
    Tandem::rebase(repeats, n_shift_map);
    
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
