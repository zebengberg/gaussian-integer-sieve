#include <iostream>
#include "../include/Moat.hpp"
#include "../include/OctantSieve.hpp"
#include "../include/OctantDonutSieve.hpp"
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




// Call this static setter method before any instances of this class are created.
void BlockMoat::setStatics(int32_t realPart, double js, bool vb) {
    BlockMoat::verbose = vb;
    BlockMoat::jumpSize = js;
    BlockMoat::dx = 100;
    BlockMoat::dy = 100;
    BlockMoat::sievingPrimesNormBound = 2 * realPart;

    // Setting nearest neighbors.
    for (int32_t u = -int32_t(jumpSize); u < jumpSize; u++) {
        for (int32_t v = -int32_t(jumpSize); v < jumpSize; v++) {
            // u and v shouldn't both be 0
            // apart from the prime 1 + i, u and v should have same parity
            // recall that c++ calculates (-3) % 2 as -1
            if ((u * u + v * v <= jumpSize * jumpSize) && (u || v) && (abs(u) % 2 == abs(v) % 2)) {
                BlockMoat::nearestNeighbors.emplace_back(u, v);
            }
        }
    }

    // Setting sieving primes.
    OctantDonutSieve d(sievingPrimesNormBound);
    d.run();
    BlockMoat::sievingPrimes = d.getBigPrimes();
}

// Constructor
BlockMoat::BlockMoat(int32_t x, int32_t y)
    // Calling BlockSieve's constructor
    : BlockSieve(x, y, dx, dy, verbose)
    , x(x)
    , y(y)
{ upperWallYPunch = dy; }

// Cannot call virtual methods of BlockSieve parent from BlockMoat constructor.
void BlockMoat::callSieve() {
    vector<gint> smallPrimes;
    for (gint g : sievingPrimes) {
        if (g.norm() <= pow((uint64_t)(x + dx - 1), 2) + pow((uint64_t)(y + dy - 1), 2)) {
            smallPrimes.push_back(g);
        }
    }
    setSmallPrimesFromReference(smallPrimes);
    setSieveArray();
    sieve();
}


// Set upperWallFlag to true if exploring from upper wall.
// If in left wall explore mode, return true if punch through the right wall
// and otherwise return false. If in upper wall explore mode, return false.
bool BlockMoat::exploreAtGint(int32_t a, int32_t b, bool upperWallFlag) {
    // Because of parity constraint in nearestNeighbors, this method will
    // never visit the ramifying prime 1 + i.
    vector<gint> toExplore;
    toExplore.emplace_back(a, b);

    do {
        gint p = toExplore.back();
        toExplore.pop_back();
        sieveArray[p.a][p.b] = false;  // indicating a visit
        for (const gint &q : nearestNeighbors) {
            gint g = p + q;

            // Checking if we've punched through
            if (upperWallFlag) {
                if (g.a >= dx) {
                    if (p.b < upperWallYPunch) {
                        upperWallYPunch = p.b;
                    }
                }
                if (g.b < 0) {  // If blocks are tall enough, this should never happen
                    cerr << "Punched through LOWER wall at: " << g.a << " " << g.b << endl;
                    cerr << "Started this exploration at: " << a << " " << b << endl;
                    // For debugging purposes
                    printSieveArray();
                    exit(1);
                }
            } else {
                if (g.a >= dx) {
                    if (verbose) {
                        cerr << "Punched through right wall at: " << g.a << " " << g.b << endl;
                        cerr << "Started this exploration at: " << a << " " << b << endl;
                        cerr << "Moving exploration block to the right..." << endl;
                    }
                    return true;
                }
            }

            // If we are not interacting with the boundary, keep exploring. The
            // check g.a < dx is somewhat redundant.
            if ((g.a >= 0) && (g.a < dx) && (g.b >= 0) && (g.b < dy) && sieveArray[g.a][g.b]) {
                    toExplore.push_back(g);
            }
        }
    }
    while (!toExplore.empty());
    return false;
}


// ALGORITHM:
// Left wall portion.
// Explore all components emanating from a gint within jumpSize of left wall
// of the current block. Use the sieve array to mark any gint that has been
// visited along the way. If we punch through the right-hand side of the block,
// immediately stop and replace current block with one to the right. If this
// punching through happens repeatedly, consider starting with a larger
// real part.
// Upper wall portion.
// Using the same component search as above, start search gints emanating from
// within jumpSize of the upper wall of the current block. We will never punch
// through the left wall because these primes have already been visited. If we
// punch through the bottom wall, immediately stop and debug. If the block
// we are searching within is much more tall than wide, this should not occur.
// Track the y-coordinates of all gints encountered within jumpSize of the right
// wall. The minimum of these y-coordinates will be the starting y-value of the
// next block to be explored.

// Return true if punch through right wall.
bool BlockMoat::exploreLeftWall() {
    for (int32_t a = 0; a < jumpSize; a++) {
        for (int32_t b = 0; b < dy; b++) {
            if (sieveArray[a][b]) {
                // Check if we've punched through block.
                if (exploreAtGint(a, b)) {
                    // Break out of this function with return.
                    return true;
                }
            }
        }
    }
    return false;
}


void BlockMoat::exploreUpperWall() {
    for (int32_t b = dy - 1; b >= dy - jumpSize; b--) {
        for (int32_t a = 0; a < dx; a++) {
            if (sieveArray[a][b]) {
                exploreAtGint(a, b, true);
            }
        }
    }
}

pair<int32_t, int32_t> BlockMoat::getNextBlock() {
    if (exploreLeftWall()) { // punched through right wall
        return {x + dx, y};
    } else {
        exploreUpperWall();
        return {x, y + upperWallYPunch};
    }
}

// Function to use BlockMoat.
void verticalMoat(int32_t realPart, double jumpSize, bool verbose) {
    BlockMoat::setStatics(realPart, jumpSize, verbose);
    int32_t x = realPart;
    int32_t y = 0;

    while(y < x) {
        // TODO: Check if this object needs to be explicitly destructed
        BlockMoat b(x, y);
        b.callSieve();
        pair<int32_t, int32_t> p = b.getNextBlock();
        x = p.first;
        y = p.second;
    }
}