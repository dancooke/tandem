all:
	g++ example.cpp tandem.cpp libdivsufsort/divsufsort.c libdivsufsort/sssort.c libdivsufsort/trsort.c libdivsufsort/utils.c -w -std=c++14 -O3 -o example