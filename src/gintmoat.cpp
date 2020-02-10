#include <iostream>
#include "../include/Moat.hpp"



int main() {
    verticalMoat(10000, 2.3);
}

int main2(int argc, const char* argv[]) {
    if (argc < 2) {
        cerr << "\n";
        cerr << "Cannot understand input. Use -h optional flag for help.\n" << endl;
        return 1;
    }

    bool origin = false;
    bool vertical = false;
    bool verbose = false;
    bool printPrimes = false;
    bool write = false;

    uint64_t x = 0;
    double jumpSize = 0;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        // Usage information.

        if ((arg == "-h") || (arg == "--help")) {
            cerr << "\n";
            cerr << "Usage: " << argv[0] << " jumpSize [x] [option1] [option2] ...\n"
                 << "Calculate the Gaussian prime moat with .\n"
                 << "    jumpSize            The jump threshold under which primes are adjacent."
                 << "    x                   The norm-bound of the search space or the \n"
                 << "                        real-part of the vertical strip to explore.\n\n"
                 << "Exploration modes:\n"
                 << "    -o, --origin        Explore the connected component of the graph starting at the\n"
                 << "                        origin. Search space holds Gaussian integers in the first\n"
                 << "                        octant. This is the default exploration mode.\n"
                 << "    --vertical          Search for Gaussian moat along a thin vertical strip starting\n"
                 << "                        at real-part x.\n\n"
                 << "Options:\n"
                 << "    -h, --help          Print this help message.\n"
                 << "    -v, --verbose       Display progress.\n"
                 << "    -p, --printprimes   Print the real and imag part of primes in the connected component\n"
                 << "                        if in origin mode."
                 << endl;
            return 1;
        }

        // Getting flags
        if ((arg == "-v") || (arg == "--verbose")) { verbose = true; }
        if ((arg == "-p") || (arg == "--printprimes")) { printPrimes = true; }
        if (arg == "--vertical") { vertical = true; }
        else { origin = true; }

        // Getting the jumpSize if it is a decimal type number.
        if (arg.front() == '.') {
            jumpSize = stod(arg);
        } else if (isdigit(arg.front())) {
            if (jumpSize == 0.0) {
                jumpSize = stod(arg);
            } else {  // Getting the parameter x.
                x = stoull(arg);
            }
        }
    }


    if (verbose) { cerr << '\n' << endl; }
    if (jumpSize == 0.0) {  // If jumpSize hasn't been parsed, abort.
        cerr << "\nCannot understand input. Use -h optional flag for help.\n" << endl;
        return 1;
    }

    if (origin) {
        if (!x) {
            x = 1000000;
        }
        OctantMoat m(x, jumpSize);
        m.exploreComponent(0, 0);
        if (printPrimes) {
            m.printCurrentComponent();
        }
        cerr << "The discovered component has size: " << m.getComponentSize() << endl;
        gint g = m.getComponentMaxElement();
        cerr << "The furthest out prime in component has coordinates: " << g.a << " " << g.b << endl;
    }
    if (vertical) {
        verticalMoat(x, jumpSize, verbose);
    }
    return 0;
}

