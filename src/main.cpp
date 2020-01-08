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
    long x = 0;
    long y = 0;
    long z = 0;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        // Usage information.
        if ((arg == "-h") || (arg == "--help")) {
            cerr << "\n";
            cerr << "usage: " << argv[0] << " x [y z] [-h] [-v] [-p] [-w] [-a] [-c] [-q] [-o] [-d] [-s]\n"
                 << "    x                   the bound on the norm of the generated primes\n"
                 << "    y                   coordinates (x, y) are lower left corner of block in segmented sieve\n"
                 << "    z                   the side length of the square block used in the segmented sieve\n\n"
                 << "options:\n"
                 << "    -h, --help          show this help message\n"
                 << "    -v, --verbose       display sieving progress\n"
                 << "    -p, --printprimes   print out the results of the sieve\n"
                 << "    -w, --write         write primes to csv file in current directory\n"
                 << "    -a, --printarray    print a text representation of the sieve array\n"
                 << "    -c, --count         count the number of generated primes and exit program\n"
                 << "    -q, --quadrant      sieve array contains all Gaussian integers in the first quadrant\n"
                 << "    -o, --octant        sieve array contains all Gaussian integers in the first octant\n"
                 << "    -d, --donut         sieve array contains Gaussian integers in first octant coprime to 2 and 5\n"
                 << "    -s, --segment       sieve array contains Gaussian integers of form a + bi with\n"
                 << "                        x <= a <= x + z and y <= b <= y + z\n"
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
// * write naive algorithm in c++
// * use primesieve with naive algorithm
// * write tests.cpp file
// * write race.cpp file
// * determine how to put c array or vector into numpy array through cython
// * use this to get results of race, possible load as np
// * rewrite readme
// * publish to pip
// * /