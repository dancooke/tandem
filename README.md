# tandem

tandem is a small C++14 library for finding tandem repeats in strings. At the moment it only finds exact repeats, but I may extend this at some point to find approximate repeats.

`main.cpp` includes a simple command line tool for DNA/RNA sequences. You can build it with CMake:

```shell
$ mkdir build && cd build
$ cmake ..
$ make
```

The command line tool reads an input stream and reports repeats:

```shell
echo NNNACGTACGTNNAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGNNNNNNNN | ./tandem
```

To find repeats in a FASTA reference genome, you can use the command line tool found in my project [bioio](https://github.com/dancooke/bioio):

```shell
./fasta human_g1k_v37.fasta Y | ./tandem
```

To report all repeats in CSV format, add `-c -a` to the command:

```shell
./fasta human_g1k_v37.fasta Y | ./tandem -c -a
```
