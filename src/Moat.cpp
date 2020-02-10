#include <iostream>
#include "../include/Moat.hpp"
#include "../include/OctantSieve.hpp"
#include "../include/OctantDonutSieve.hpp"
#include "../include/BlockSieve.hpp"
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

vector<vector<gint>> OctantMoat::getAllComponents() {
    return allComponents;
}




VerticalMoat::VerticalMoat(uint32_t realPart, double jumpSize, bool verbose)
    : realPart(realPart)
    , jumpSize(jumpSize)
    , verbose(verbose)
{
    normBound = 2 * realPart;
    OctantDonutSieve d(normBound);
    d.run();
    sievingPrimes = d.getBigPrimes();  // to be passed into all instances of BlockSieve
    x = realPart - (realPart % 10);  // must be a multiple of 10
    y = 0;
    dx = 100;
    dy = 20;
    setNearestNeighbors();
}

void VerticalMoat::setNearestNeighbors() {
    for (int32_t u = -int32_t(jumpSize); u < jumpSize; u++) {
        for (int32_t v = -int32_t(jumpSize); v < jumpSize; v++) {
            // u and v shouldn't both be 0
            // apart from the prime 1 + i, u and v should have same parity
            // recall that c++ calculates (-3) % 2 as -1
            if ((u * u + v * v <= jumpSize * jumpSize) && (u || v) && (abs(u) % 2 == abs(v) % 2)) {
                VerticalMoat::nearestNeighbors.emplace_back(u, v);
            }
        }
    }
}

void VerticalMoat::setSieveArray(){
    BlockSieve b(x, y, dx, dy);
    vector<gint> smallPrimes;
    for (gint g : sievingPrimes) {
        if (g.norm() <= pow((uint64_t)(x + dx - 1), 2) + pow((uint64_t)(y + dy - 1), 2)) {
            smallPrimes.push_back(g);
        }
    }
    b.setSmallPrimesFromReference(smallPrimes);
    b.setSieveArray();
    b.sieve();
    sieveArray = b.getSieveArray();
}

void VerticalMoat::printSieveArray() {
    for (int32_t b = dy - 1; b >= 0; b--) {
        string row;
        for (int32_t a = 0; a < dx; a++) {
            row += sieveArray[a][b] ? '*' : '-';
        }
        cerr << row << endl;
    }
}

int VerticalMoat::exploreAtGint(int32_t a, int32_t b) {
    // The gint starting_g must be prime; ignoring 1 + i
    gint starting_g(a, b);
    // reset current component
    currentComponent.clear();
    vector<gint> toExplore;
    toExplore.push_back(starting_g);
    currentComponent.push_back(starting_g);
    sieveArray[starting_g.a][starting_g.b] = false;

    uint32_t count = 0;
    while (!toExplore.empty()) {
        gint p = toExplore.back();
        toExplore.pop_back();
        for (const gint &q : nearestNeighbors) {
            gint g = p + q;
            // Checking if inside block and prime
            if (g.a >= dx) {
                if (verbose) {
                    cerr << "We have punched through the block at: " << g.a << " " << g.b << endl;
                    cerr << "We started this exploration at: " << a << " " << b << endl;
                }
                return 1;
            } else if ((g.a >= 0) && (g.b >= 0) && (g.b < dy) && sieveArray[g.a][g.b]) {
                currentComponent.push_back(g);
                sieveArray[g.a][g.b] = false;  // indicating that g has been visited
                toExplore.push_back(g);
            }
        }
    }
    return 0;
}

void VerticalMoat::exploreLeftWall() {
    for (int32_t a = 0; a < jumpSize; a++) {
        for (int32_t b = 0; b < dy; b++) {
            if (sieveArray[a][b]) {
                if (exploreAtGint(a, b)) {  // check if we've punched through block
                    // force a break out of double for loop
                    b = dy;
                    a = jumpSize;
                }
            }
        }
    }
}

void VerticalMoat::exploreUpperWall() {
    for (int32_t b = dy - 1; b >= dy - jumpSize; b--) {
        for (int32_t a = 0; a < dx; a++) {
            if (sieveArray[a][b]) {
                exploreAtGint(a, b);
            }
        }
    }
}

