#include <iostream>
#include "Moat.hpp"


int main() {
    uint64_t s = getCount();
    cout << "\n\nThe size of the component: " << s << endl;
}




int main2(int argc, const char* argv[]) {
    if (argc < 3) {
        cerr << "\n";
        cerr << "Not enough parameters passed. Use -h optional flag for help.\n" << endl;
        return 1;
    }

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
            cerr << "Usage: " << argv[0] << " x jumpSize [option1] [option2] ...\n"
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

        // Parsing for jumpSize and x.
        if (isdigit(arg.front())) {
            if (!x) {
                x = stoull(arg);
            } else {
                jumpSize = stod(arg);
            }
        }
    }


    if (verbose) { cerr << '\n' << endl; }
    if ((jumpSize == 0.0) || !x) {  // If jumpSize hasn't been parsed, abort.
        cerr << "\nCannot understand input. Use -h optional flag for help.\n" << endl;
        return 1;
    }


    if (vertical) {
        if (verbose) {
            cerr << "Searching for moat in vertical strip..." << endl;
        }
        verticalMoat(x, jumpSize, verbose);
    } else {
        if (verbose) {
            cerr << "Searching for moat starting at origin..." << endl;
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

    return 0;
}

