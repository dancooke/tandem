# tandem

Tandem is a small C++ library for finding tandem repeats in strings. At the moment it only finds exact repeats, but I may extend this at some point to find approximate repeats.

example.cpp includes a basic example for DNA/RNA sequences. To compile this example just `make` in the tandem directory. You can run the example like

    echo "NNNACGTACGTNNAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGNNNNNNNN" | ./example
