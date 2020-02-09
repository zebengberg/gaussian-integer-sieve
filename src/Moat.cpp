#include <iostream>
#include "../include/Moat.hpp"
#include "../include/OctantSieve.hpp"
using namespace std;



// Public methods in OctantMoat class.
OctantMoat::OctantMoat(uint64_t normBound, double jumpSize, bool verbose)
    : normBound(normBound)
    , jumpSize(jumpSize)
    , verbose(verbose)
{
    // using tolerance with jumpSize
    double tolerance = pow(10, -3);
    this->jumpSize += tolerance;
    OctantSieve o(normBound, verbose);
    o.run();
    sieveArray = o.getSieveArray();
    setNearestNeighbors();
}


void OctantMoat::setNearestNeighbors() {
    for (int32_t u = -int32_t(jumpSize); u < jumpSize; u++) {
        for (int32_t v = -int32_t(jumpSize); v < jumpSize; v++) {
            // u and v shouldn't both be 0
            // apart from the prime 1 + i, u and v should have same parity
            // recall that c++ calculates (-3) % 2 as -1
            if ((u * u + v * v <= jumpSize * jumpSize) && (u || v) && (abs(u) % 2 == abs(v) % 2)) {
                nearestNeighbors.emplace_back(u, v);
            }
        }
    }
}

// Depth first search exploring the connected component of starting_g
void OctantMoat::exploreComponent(int32_t a, int32_t b) {
    gint starting_g(a, b);
    // reset current component
    currentComponent.clear();
    vector<gint> toExplore;


    // The gint starting_g must either be zero or a prime.
    if (starting_g.a || starting_g.b) {  // prime
        toExplore.push_back(starting_g);
        currentComponent.push_back(starting_g);
        sieveArray[starting_g.a][starting_g.b] = false;
    } else {  // zero
        if (jumpSize > sqrt(2)) {
            currentComponent.emplace_back(1, 1);
            currentComponent.emplace_back(2, 1);
            sieveArray[1][1] = false;
            sieveArray[2][1] = false;
            toExplore.emplace_back(2, 1);
        }
    }

    uint32_t count = 0;
    while (!toExplore.empty()) {
        gint p = toExplore.back();
        toExplore.pop_back();
        for (const gint &q : nearestNeighbors) {
            gint g = p + q;
            if (g.norm() <= normBound) {
                // Checking if inside first octant and prime
                if ((g.a >= 0) && (g.b >= 0) && (g.b <= g.a) && sieveArray[g.a][g.b]) {
                    currentComponent.push_back(g);
                    sieveArray[g.a][g.b] = false;  // indicating that g has been visited
                    toExplore.push_back(g);
                }
            } else {
                // Checking if we have punched through without encountering a moat when starting at 0, 0
                if (!starting_g.a && !starting_g.b) {
                    cerr << "\nTraversed outside of the norm bound!" << endl;
                    cerr << "Failed to find a moat of size " << jumpSize << endl;
                    exit(1);
                }
            }
        }
        count++;
        if (verbose) {
            if (count % 1000 == 0) {
                cerr << '.';
            }
            if (count % 100000 == 0) {
                cerr << endl;
            }
        }
    }
    cout << endl;
}

uint32_t OctantMoat::getComponentSize() {
    return currentComponent.size();
}

gint OctantMoat::getComponentMaxElement() {
    return *max_element(currentComponent.begin(), currentComponent.end());
};

vector<gint> OctantMoat::getCurrentComponent() {
    return currentComponent;
}

void OctantMoat::printCurrentComponent() {
    for (gint g : currentComponent) {
        cout << g.a << " " << g.b << endl;
    }
}

// The sieve array was modified in explore(), and this method gets primes still marked true.
vector<gint> OctantMoat::getUnexplored() {
    vector<gint> unexplored;
    for (uint32_t u = 0; u < sieveArray.size(); u++) {
        for (uint32_t v = 0; v < sieveArray[u].size(); v++) {
            if (sieveArray[u][v]) {
                unexplored.emplace_back(u, v);
            }
        }
    }
    return unexplored;
}

void OctantMoat::exploreAllComponents() {
    for (uint32_t u = 0; u < sieveArray.size(); u++) {
        for (uint32_t v = 0; v < sieveArray[u].size(); v++) {
            if (sieveArray[u][v]) {
                exploreComponent(u, v);
                allComponents.push_back(currentComponent);
            }
        }
    }
}

// Python will convert this to a list of pointers to arrays.
vector<vector<gint>> OctantMoat::getAllComponents() {
    return allComponents;
}


int main(int argc, const char* argv[]) {
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
        if (jumpSize == 0.0) {
            if (isdigit(arg.front()) || (arg.front() == '.')) {
                jumpSize = stod(arg);
            }
        } else {  // Getting the parameter x.
            if (isdigit(arg.front())) {
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
        cerr << "\nExploring the connected component starting at the origin.\n" << endl;
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
        cerr << "\nNot yet implemented!" << endl;
    }
    return 0;
}

