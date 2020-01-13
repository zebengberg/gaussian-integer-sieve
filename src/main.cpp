#include <iostream>
#include "../include/QuadrantSieve.hpp"
#include "../include/OctantSieve.hpp"
#include "../include/DonutSieve.hpp"
using namespace std;

int main(int argc, const char* argv[]) {
    if (argc < 2) {
        cerr << "\n";
        cerr << "Cannot understand input. Use -h optional flag for help.\n" << endl;
        return 1;
    }

    bool verbose = false;
    bool printPrimes = false;
    bool printArray = false;
    bool write = false;
    bool count = false;
    string sieveType = "donut";  // put in default sieve method here
    uint64_t x = 0;
    uint32_t y = 0;
    uint32_t z = 0;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        // Usage information.
        if ((arg == "-h") || (arg == "--help")) {
            cerr << "\n";
            cerr << "Usage: " << argv[0] << " x [y z] [option1] [option2] ...\n"
                 << "Generate Gaussian primes with norm up to x using sieving methods.\n"
                 << "    x                   The norm-bound of the generated primes\n"
                 << "    y                   Coordinates (x, y) of SW-corner of block in segmented sieve.\n"
                 << "    z                   Side length block in segmented sieve.\n\n"
                 << "Options:\n"
                 << "    -h, --help          Print this help message.\n"
                 << "    -v, --verbose       Display sieving progress.\n"
                 << "    -p, --printprimes   Print the real and imag part of primes found by the sieve.\n"
                 << "    -w, --write         Write primes to csv file in working directory.\n"
                 << "    -a, --printarray    Print a text representation of the sieve array.\n"
                 << "    -c, --count         Count the number of generated primes and exit program.\n"
                 << "    -q, --quadrant      Sieve array consists of Gaussian integers in the first quadrant.\n"
                 << "    -o, --octant        Sieve array consists of Gaussian integers in the first octant.\n"
                 << "    -d, --donut         Sieve array consists of Gaussian integers in first octant\n"
                 << "                        coprime to 2 and 5.\n"
                 << "    -s, --segmented     Sieve array consists of Gaussian integers of form a + bi with\n"
                 << "                        x <= a < x + z and y <= b < y + z.\n"
                 << std::endl;
            return 1;
        }

        // Getting flags
        if ((arg == "-v") || (arg == "--verbose")) { verbose = true; }
        if ((arg == "-p") || (arg == "--printprimes")) { printPrimes = true; }
        if ((arg == "-w") || (arg == "--write")) { write = true; }
        if ((arg == "-a") || (arg == "--printarray")) { printArray = true; }
        if ((arg == "-c") || (arg == "--count")) { count = true; }

        if ((arg == "-q") || (arg == "--quadrant")) { sieveType = "quadrant"; }
        else if ((arg == "-o") || (arg == "--octant")) { sieveType = "octant"; }
        else if ((arg == "-s") || (arg == "--segmented")) { sieveType = "segmented"; }

        // Use first char of string to determine if string is a number.
        if (isdigit(arg.front())) {
            if (y) {
                z = stol(arg);
            } else if (x) {
                y = stol(arg);
            } else {
                x = stol(arg);
            }
        }
    }
    // If x hasn't been parsed, abort.
    if (!x) {
        cerr << "\n";
        cerr << "Cannot understand input. Use -h optional flag for help.\n" << endl;
        return 1;
    }

    if (verbose) { cerr << '\n' << endl; }
    // Very boilerplate, but sieve objects are distinct.
    if (sieveType == "donut") {
        DonutSieve s(x, verbose);
        s.run();
        if (printArray) { s.printSieveArray(); }
        if (count) {
            s.getCountBigPrimes();
            return 0; // early exit for count
        }
        s.setBigPrimes();
        s.sortBigPrimes();
        if (write) {
            s.writeBigPrimesToFile();
        }
        // Default behavior if no useful options passed in.
        if (printPrimes || ((!printPrimes) && (!printArray) && (!write))) {
            s.printBigPrimes();
        }
    } else if (sieveType == "quadrant") {
        QuadrantSieve s(x, verbose);
        s.run();
        if (printArray) { s.printSieveArray(); }
        if (count) {
            s.getCountBigPrimes();
            return 0;  // early exit for count
        }
        s.setBigPrimes();
        s.sortBigPrimes();
        if (write) {
            s.writeBigPrimesToFile();
            return 0;
        }
        // Default behavior if no useful options passed in.
        if (printPrimes || ((!printPrimes) && (!printArray) && (!write))) {
            s.printBigPrimes();
        }
    } else if (sieveType == "octant") {
        OctantSieve s(x, verbose);
        s.run();
        if (printArray) { s.printSieveArray(); }
        if (count) {
            s.getCountBigPrimes();
            return 0;  // early exit for count
        } else {
            s.setBigPrimes();
            s.sortBigPrimes();
        }
        if (write) {
            s.writeBigPrimesToFile();
        }
        // Default behavior if no useful options passed in.
        if (printPrimes || ((!printPrimes) && (!printArray) && (!write))) {
            s.printBigPrimes();
        }
    } else {
        SegmentedSieve s(x, y, z, verbose);
        s.run();
        if (printArray) { s.printSieveArray(); }
        if (count) {
            s.getCountBigPrimes();
            return 0;  // early exit for count
        }
        s.setBigPrimes();
        s.sortBigPrimes();
        if (write) {
            s.writeBigPrimesToFile();
        }
        // Default behavior if no useful options passed in.
        if (printPrimes || ((!printPrimes) && (!printArray) && (!write))) {
            s.printBigPrimes();
        }
    }
    return 0;
}



///* TODO
// * write tests.cpp file
// * write race.cpp file
// * determine how to put c array or vector into numpy array through cython
// * use this to get results of race, possible load as np
// * publish to pip
// * /