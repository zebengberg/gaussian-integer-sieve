#include <iostream>
#include "../include/Moat.hpp"
#include "../include/OctantSieve.hpp"
using namespace std;



// Public methods in OctantMoat class.
OctantMoat::OctantMoat(double jumpSize, bool verbose)
    : jumpSize(jumpSize)
    , verbose(verbose)
{
    // using tolerance with jumpSize
    double tolerance = pow(10, -3);
    this->jumpSize += tolerance;

    // setting normBound according to existing tables in Tsuchimura paper
    if (jumpSize < 2.1) {  // jumpsize <= 2
        normBound = 3000;
    } else if (jumpSize < 3) {  // jumpsize = sqrt(8)
        normBound = 10000;
    } else if (jumpSize < 4) {  // jumpsize = sqrt(10)
        normBound = 1100000;
    } else if (jumpSize < 4.2) {  // jumpsize = 4
        normBound = 20000000;
    } else if (jumpSize < 4.4) {  // jumpsize = sqrt(18)
        normBound = 116000000;
    } else if (jumpSize < 5) {  // jumpsize = sqrt(20)
        normBound = 17900000000;
    } else {
        cerr << "Jump size is too large for this method!" << endl;
        cerr << "Instead call the segmented moat" << endl;
        exit(1);
    }


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
            if (u * u + v * v <= jumpSize * jumpSize && (u || v) && abs(u) % 2 == abs(v) % 2) {
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
        if (verbose) {
            cerr << "\nExploring the connected component starting at the origin.\n" << endl;
        }
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
            if (count % 10000 == 0) {
                cerr << '.';
            }
            if (count % (80 * 10000) == 0) {
                cerr << endl;
            }
        }
    }
    cerr << endl;
}

uint32_t OctantMoat::getComponentSize() {
    return currentComponent.size();
}

gint OctantMoat::getComponentMaxElement() {
    return *max_element(currentComponent.begin(), currentComponent.end());
}

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

vector<vector<gint>> OctantMoat::getAllComponents() {
    return allComponents;
}



