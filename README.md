# tandem

Tandem is a small C++ library for finding tandem repeats in strings. At the moment it only finds exact repeats, but I may extend this at some point to find approximate repeats.

example.cpp includes a simple example for DNA/RNA sequences. To compile the example just `make` in the tandem directory. You can run it like

    echo "NNNACGTACGTNNAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGNNNNNNNN" | ./example

