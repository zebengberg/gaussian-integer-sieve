#include <iostream>
#include <random>
#include "../include/OctantSieve.hpp"
#include "../include/OctantDonutSieve.hpp"
#include "../include/BlockSieve.hpp"
#include "../include/BlockDonutSieve.hpp"
#include "../include/SectorSieve.hpp"
using namespace std;

// Delete after debugging
int main() {
    uint64_t x = pow(2, 25);
    random_device rd;
    uniform_real_distribution<double> distReal(0.0, M_PI_4 - 0.01);
    while (true) {
        double alpha = distReal(rd);
        double beta = distReal(rd);
        if (beta < alpha) {
            double temp = beta;
            beta = alpha;
            alpha = temp;
        }
        cout << alpha << " " << beta << endl;

        SectorSieve s(x, alpha, beta, false);
        s.run();
        vector<gint> sP = s.getBigPrimes(false);  // not sorting yet

        OctantDonutSieve d(x, false);
        d.run();
        vector<gint> dP = d.getBigPrimes();

        // Removing gints outside of sector
        auto res = remove_if(dP.begin(), dP.end(), [&](gint g) { return (g.arg() >= beta) || (g.arg() < alpha); });
        dP.erase(res, dP.end());

        // Sorting generated primes and checking if two lists are equal.
        sort(sP.begin(), sP.end());
        sort(dP.begin(), dP.end());
        assert(sP == dP);
    }
}


int main2(int argc, const char* argv[]) {
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
    bool donut = false;
    bool octant = false;
    bool block = false;
    bool sector = false;

    uint64_t x = 0;
    uint64_t y = 0;
    uint32_t dx = 0;
    uint32_t dy = 0;
    double alpha = -1.0;
    double beta = -1.0;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        // Usage information.
        if ((arg == "-h") || (arg == "--help")) {
            cerr << "\n";
            cerr << "Usage: " << argv[0] << " x [y dx dy alpha beta] [option1] [option2] ...\n"
                 << "Generate Gaussian primes with norm up to x using sieving methods.\n"
                 << "    x                   The norm-bound of the generated primes\n"
                 << "    y                   Coordinates (x, y) of SW-corner of array in block sieve.\n"
                 << "    dx                  Horizontal side length in block sieve.\n"
                 << "    dy                  Vertical side length in block sieve.\n"
                 << "    alpha               Start angle in sector sieve.\n"
                 << "    beta                Final angle in sector sieve.\n\n"
                 << "Options:\n"
                 << "    -h, --help          Print this help message.\n"
                 << "    -v, --verbose       Display sieving progress.\n"
                 << "    -p, --printprimes   Print the real and imag part of primes found by the sieve.\n"
                 << "    -w, --write         Write primes to csv file in current directory.\n"
                 << "    -a, --printarray    Print a text representation of the sieve array.\n"
                 << "    -c, --count         Count the number of generated primes and exit program.\n\n"
                 << "Optional sieve types:\n"
                 << "    -o, --octant        Sieve array indexed by Gaussian integers in the first octant.\n"
                 << "                        This is the default sieve method called.\n"
                 << "    -s, --sector        Sieve array indexed by Gaussian integers in the sector\n"
                 << "                        with start angle alpha and final angle beta.\n"
                 << "    -b, --block         Sieve array indexed by Gaussian integers in the rectangle\n"
                 << "                        defined by x <= real < x + dx and y <= imag < y + dy.\n"
                 << "    -d, --donut         If a donut version of the sieve array exists, use it. In the\n"
                 << "                        donut sieve, the sieve array consists of Gaussian integers\n"
                 << "                        coprime to 2 and 5. This optional can be used with --octant\n"
                 << "                        and --block, and is often significantly faster.\n"
                 << std::endl;
            return 1;
        }

        // Getting flags
        if ((arg == "-v") || (arg == "--verbose")) { verbose = true; }
        if ((arg == "-p") || (arg == "--printprimes")) { printPrimes = true; }
        if ((arg == "-w") || (arg == "--write")) { write = true; }
        if ((arg == "-a") || (arg == "--printarray")) { printArray = true; }
        if ((arg == "-c") || (arg == "--count")) { count = true; }
        if ((arg == "-d") || (arg == "--donut")) { donut = true; }
        if ((arg == "-o") || (arg == "--octant")) { octant = true; }
        if ((arg == "-s") || (arg == "--sector")) { sector = true; }
        if ((arg == "-b") || (arg == "--block")) { block = true; }

        // Getting the input if it is a decimal type number.
        if ((arg.front() == '0') || (arg.front() == '.')) {
            if (!(alpha == -1.0)) {
                beta = stod(arg);
            } else {
                alpha = stod(arg);
            }
        }

        // Getting the input if it is an integer.
        if (isdigit(arg.front())) {
            if (dx) {
                dy = stoul(arg);
            } else if (y) {
                dx = stoul(arg);
            } else if (x) {
                y = stoull(arg);
            } else {
                x = stoull(arg);
            }
        }
    }

    // Getting sieve type.
    string sieveType = "octantDonut";  // put in default here
    if (!x) {  // If x hasn't been parsed, abort.
        cerr << "\nCannot understand input. Use -h optional flag for help.\n" << endl;
        return 1;
    } else {
        if (sector) {
            sieveType = "sector";
        } else if (block && donut) {
            sieveType = "blockDonut";
        } else if (block) {
            sieveType = "block";
        } else if (octant && !donut) {
            sieveType = "octant";
        }
    }

    if (verbose) { cerr << '\n' << endl; }
    // Very boilerplate, but sieve objects are distinct.
    if (sieveType == "octantDonut") {
        cerr << "\nCalling the Octant Donut Sieve.\n" << endl;
        OctantDonutSieve s(x, verbose);
        s.run();
        if (printArray) { s.printSieveArray(); }
        if (count) {
            cout << s.getCountBigPrimes() << endl;
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
    } else if (sieveType == "octant") {
        cerr << "\nCalling the Octant Sieve.\n" << endl;
        OctantSieve s(x, verbose);
        s.run();
        if (printArray) { s.printSieveArray(); }
        if (count) {
            cout << s.getCountBigPrimes() << endl;
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
    } else if (sieveType == "sector"){
        if ((alpha == -1.0) || (beta == -1.0)) {
            cerr << "Provide angle values to use sector sieve.\n"
                 << "Use -h optional flag for help.\n" << endl;
            return 1;
        }
        cerr << "\nCalling the Sector Sieve.\n" << endl;
        SectorSieve s(x, alpha, beta, verbose);
        s.run();
        if (printArray) { s.printSieveArray(); }
        if (count) {
            cout << s.getCountBigPrimes() << endl;
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
    } else if (sieveType == "blockDonut") {
        if (!y || !dx || !dy) {
            cerr << "Provide coordinates x, y, dx, and dy to use block sieve.\n"
                 << "Use -h optional flag for help.\n" << endl;
            return 1;
        }
        cerr << "\nCalling the Block Donut Sieve.\n" << endl;
        BlockDonutSieve s(x, y, dx, dy);
        s.run();
        if (printArray) { s.printSieveArray(); }
        if (count) {
            cout << s.getCountBigPrimes() << endl;
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
    else if (sieveType == "block") {
        if (!y || !dx || !dy) {
            cerr << "Provide coordinates x, y, dx, and dy to use block sieve.\n"
                 << "Use -h optional flag for help.\n" << endl;
            return 1;
        }
        cerr << "\nCalling the Block Sieve.\n" << endl;
        BlockSieve s(x, y, dx, dy);
        s.run();
        if (printArray) { s.printSieveArray(); }
        if (count) {
            cout << s.getCountBigPrimes() << endl;
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
// * write tests.cpp file and benchmarks file
// * write race.cpp file
// * use this to get results of race, possible load as np
// * publish to pip
// rewrite readme
// * /