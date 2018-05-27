#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <map>

#include "tandem.hpp"

/**
    Simple example showing how to find all exact tandem repeats in a DNA sequence using
    the Tandem API.
 */

bool cmd_option_exists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

std::string get_cmd_option(char** begin, char** end, const std::string& option)
{
    auto it = std::find(begin, end, option);
    if (it != end && ++it != end) return *it;
    return {};
}

template <typename T>
std::string to_string(const T val, const unsigned precision = 2)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << val;
    return out.str();
}

namespace detail {

struct CapitaliseBase
{
    auto operator()(const char base) const noexcept
    {
        switch (base) {
            case 'a': return 'A';
            case 'c': return 'C';
            case 'g': return 'G';
            case 't': return 'T';
            case 'u': return 'U';
            case 'n': return 'N';
            default : return base;
        }
    }
};

} // namespace detail

template <typename SequenceType>
void capitalise(SequenceType& sequence)
{
    std::transform(std::begin(sequence), std::end(sequence), std::begin(sequence),
                   detail::CapitaliseBase {});
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

std::vector<tandem::Repeat> get_repeats(std::string sequence)
{
    if (sequence.empty()) return {};
    std::clog << "Looking for repeats in " << sequence.size() << "bp" << std::endl;
    if (sequence.back() != '$') {
        assert(std::find(std::cbegin(sequence), std::cend(sequence), '$') == std::cend(sequence));
        sequence.reserve(sequence.size() + 1);
        sequence.push_back('$'); // last character must not be in alphabet
    }
    auto n_shift_map = tandem::collapse(sequence, 'N'); // makes search a lot faster
    auto repeats = tandem::extract_exact_tandem_repeats(sequence);
    tandem::rebase(repeats, n_shift_map);
    return repeats;
}

struct PeriodGreater
{
    bool operator()(const tandem::Repeat& lhs, const tandem::Repeat& rhs) const noexcept
    {
        return lhs.period > rhs.period;
    }
    bool operator()(const std::size_t& lhs, const tandem::Repeat& rhs) const noexcept
    {
        return lhs > rhs.period;
    }
    bool operator()(const tandem::Repeat& lhs, const std::size_t& rhs) const noexcept
    {
        return lhs.period > rhs;
    }
};

struct LengthGreater
{
    bool operator()(const tandem::Repeat& lhs, const tandem::Repeat& rhs) const noexcept
    {
        return lhs.length > rhs.length;
    }
    bool operator()(const std::size_t& lhs, const tandem::Repeat& rhs) const noexcept
    {
        return lhs > rhs.length;
    }
    bool operator()(const tandem::Repeat& lhs, const std::size_t& rhs) const noexcept
    {
        return lhs.length > rhs;
    }
};

struct PrintableRepeat
{
    const std::string& sequence;
    const tandem::Repeat& repeat;
};

auto pos(const PrintableRepeat& pr) noexcept { return pr.repeat.pos; }
auto period(const PrintableRepeat& pr) noexcept { return pr.repeat.period; }
auto length(const PrintableRepeat& pr) noexcept { return pr.repeat.length; }
auto end(const PrintableRepeat& pr) noexcept {return pos(pr) + length(pr); }

auto num_periods(const tandem::Repeat& repeat) noexcept
{
    return static_cast<double>(repeat.length) / repeat.period;
}

auto num_periods(const PrintableRepeat& pr) noexcept
{
    return num_periods(pr.repeat);
}

std::ostream& operator<<(std::ostream& os, const PrintableRepeat& repeat)
{
    os << "@ " << pos(repeat) << '-' << end(repeat) << ' ';
    const auto first_base_itr = std::next(std::cbegin(repeat.sequence), pos(repeat));
    std::copy(first_base_itr, std::next(first_base_itr, period(repeat)), std::ostreambuf_iterator<char> {os});
    os << " ; " << period(repeat) << " x " << to_string(num_periods(repeat)) << " = " << length(repeat);
    return os;
}

using RepeatIterator = std::vector<tandem::Repeat>::iterator;

void print(RepeatIterator first, RepeatIterator last, const std::string& sequence)
{
    std::transform(first, last, std::ostream_iterator<PrintableRepeat> {std::cout, "\n"},
                   [&] (const auto& repeat) -> PrintableRepeat { return {sequence, repeat}; });
}

void report_biggest_periods(std::vector<tandem::Repeat>& repeats, const std::string& sequence,
                            const std::size_t max, const int min_size = -1)
{
    auto n = std::min(max, repeats.size());
    auto nth = std::next(std::begin(repeats), n);
    std::partial_sort(std::begin(repeats), nth, std::end(repeats), PeriodGreater {});
    if (min_size > 0) {
        nth = std::lower_bound(std::begin(repeats), nth, static_cast<std::size_t>(min_size), PeriodGreater {});
        n = std::distance(std::begin(repeats), nth);
    }
    std::clog << "Printing the " << n << " repeats with the largest periods" << std::endl;
    print(std::begin(repeats), nth, sequence);
}

void report_biggest_lengths(std::vector<tandem::Repeat>& repeats, const std::string& sequence,
                            const std::size_t max, const int min_size = -1)
{
    auto n = std::min(max, repeats.size());
    auto nth = std::next(std::begin(repeats), n);
    std::partial_sort(std::begin(repeats), nth, std::end(repeats), LengthGreater {});
    if (min_size > 0) {
        nth = std::lower_bound(std::begin(repeats), nth, static_cast<std::size_t>(min_size), LengthGreater {});
        n = std::distance(std::begin(repeats), nth);
    }
    std::clog << "Printing the " << n << " repeats with the largest lengths" << std::endl;
    print(std::begin(repeats), nth, sequence);
}

int main(int argc, char** argv)
{
    std::string sequence;
    std::getline(std::cin, sequence);
    capitalise(sequence);
    if (!is_dna(sequence) && !is_rna(sequence)) {
        std::clog << "this example is only for DNA or RNA sequences" << std::endl;
        exit(0);
    }
    auto repeats = get_repeats(sequence);
    std::size_t max_report {10};
    if (cmd_option_exists(argv, argv + argc, "-n")) {
        const std::string user_max_report {get_cmd_option(argv, argv + argc, "-n")};
        max_report = std::stoull(user_max_report);
    }
    std::clog << "Found " << repeats.size() << " tandem repeats!" << std::endl;
    int min_period {-1};
    if (cmd_option_exists(argv, argv + argc, "-p")) {
        const std::string user_min_period {get_cmd_option(argv, argv + argc, "-p")};
        min_period = std::stoi(user_min_period);
    }
    //report_biggest_periods(repeats, sequence, max_report, min_period);
    int min_length {-1};
    if (cmd_option_exists(argv, argv + argc, "-l")) {
        const std::string user_min_length {get_cmd_option(argv, argv + argc, "-l")};
        min_length = std::stoi(user_min_length);
    }
    report_biggest_lengths(repeats, sequence, max_report, min_length);
}
