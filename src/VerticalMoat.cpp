// ALGORITHM:
// Left wall portion.
// Explore all components emanating from a gint within jumpSize of left wall
// of the current block. Use the sieve array to mark any gint that has been
// visited along the way. If we punch through the right-hand side of the block,
// immediately stop and replace current block with one to the right. If this
// punching through happens repeatedly, consider starting with a larger
// real part.
//
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



#include <iostream>
#include "../include/Moat.hpp"
#include "../include/OctantDonutSieve.hpp"

// Need to first declare static member variables here in source.
bool VerticalMoat::verbose;
double VerticalMoat::jumpSize;
int32_t VerticalMoat::realPart;
int32_t VerticalMoat::blockSize;
int32_t VerticalMoat::dx;
int32_t VerticalMoat::dy;
uint64_t VerticalMoat::sievingPrimesNormBound;
vector<gint> VerticalMoat::sievingPrimes;
vector<gint> VerticalMoat::nearestNeighbors;

// Call this static setter method before any instances of this class are created.
void VerticalMoat::setStatics(int32_t rp, double js, bool vb) {
    if (vb) {
        cerr << "Setting static variables..." << endl;
    }
    realPart = rp;
    verbose = vb;
    jumpSize = js;
    // arbitrary initial values for dx and dy; these will be updated below.
    blockSize = pow(10, 7);
    dx = 1000;
    dy = blockSize / dx;
    // Bound on the norm of pre-computed primes. The factor 1.2 gives some
    // wiggle room in case there are many moves to the right.
    sievingPrimesNormBound = uint64_t(1.2 * (sqrt(2) * realPart + dx * dy));

    // Setting nearest neighbors.
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

    if (verbose) {
        cerr << "Precomputing sieving primes." << endl;
    }
    OctantDonutSieve d(sievingPrimesNormBound, false);  // not letting this be verbose
    d.run();
    sievingPrimes = d.getBigPrimes();
}

// Constructor
VerticalMoat::VerticalMoat(int32_t x, int32_t y)
// Calling BlockSieve's constructor
        : BlockSieve(x, y, dx, dy, false) // not letting this be verbose
        , x(x)
        , y(y)
{
    upperWallYPunch = dy;
    countVisited = 0;
    farthestRight = 0;
    if (verbose) {
        cerr << "Working within block having lower left corner at: " << x << " " << y << endl;
        cerr << "The block has dimensions: " << dx << " by " << dy << endl;
    }
}

// Cannot call virtual methods of parent from VerticalMoat constructor.
void VerticalMoat::callSieve() {
    // Checking to make sure there are enough primes within sievingPrimes
    gint last_g = sievingPrimes.back();
    if (last_g.norm() < isqrt(maxNorm)) {
        cerr << "Not enough pre-computed primes in static variable sievingPrimes..." << endl;
        exit(1);
    }

    vector<gint> smallPrimes;
    for (gint g : sievingPrimes) {
        if (g.norm() <= maxNorm) {
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
bool VerticalMoat::exploreAtGint(int32_t a, int32_t b, bool upperWallFlag) {
    // Because of parity constraint in nearestNeighbors, this method will
    // never visit the ramifying prime 1 + i.
    vector<gint> toExplore;
    toExplore.emplace_back(a, b);

    do {
        gint p = toExplore.back();
        toExplore.pop_back();
        sieveArray[p.a][p.b] = false;  // indicating a visit
        countVisited++;
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
                        cerr << "Moving exploration block to the right...\n" << endl;
                    }
                    return true;
                }
            }

            // If we are not interacting with the boundary, keep exploring. The
            // check g.a < dx is somewhat redundant.
            if (g.a >= 0 && g.a < dx && g.b >= 0 && g.b < dy && sieveArray[g.a][g.b]) {
                toExplore.push_back(g);
                if ((!upperWallFlag) && (g.a > farthestRight)) {
                    farthestRight = g.a;
                }
            }
        }
    }
    while (!toExplore.empty());
    return false;
}


bool VerticalMoat::exploreLeftWall() {
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


void VerticalMoat::exploreUpperWall() {
    for (int32_t b = dy - 1; b >= dy - jumpSize; b--) {
        for (int32_t a = 0; a < dx; a++) {
            if (sieveArray[a][b]) {
                exploreAtGint(a, b, true);
            }
        }
    }
}

pair<int32_t, int32_t> VerticalMoat::getNextBlock() {
    if (exploreLeftWall()) { // punched through right wall
        // double dx
        dx *= 2;
        dx = min(dx, 1000);  // clip at 1000
        dy = blockSize / dx;
        return {x + dx, y};
    } else {
        exploreUpperWall();
        // Recalibrate dx and dy passed on farthestRight
        if (farthestRight < dx / 2) {
            dx = 2 * farthestRight;
            dy = blockSize / dx;
        }
        if (verbose) {
            cerr << "Largest real-part reached starting from left hand wall: " << farthestRight << endl;
            cerr << "Number of visited primes: " << countVisited << "\n" << endl;
        }
        return {x, y + upperWallYPunch};
    }
}

// Static method to control instances of VerticalMoat.
void VerticalMoat::findVerticalMoat() {
    int32_t x = realPart;
    int32_t y = 0;
    int32_t consecutiveStepsRight = 0;

    while (y < x) {
        VerticalMoat b(x, y);
        b.callSieve();
        pair<int32_t, int32_t> p = b.getNextBlock();
        if (p.first != x) {
            consecutiveStepsRight++;
        } else {
            consecutiveStepsRight = 0;
        }
        if (consecutiveStepsRight > 10) {
            cerr << "Stepped right ten times in a row!" << endl;
            cerr << "Choose a larger value for realPart" << endl;
            exit(1);
        }
        x = p.first;
        y = p.second;
    }
    cerr << "Gaussian moat present from real-axis to boundary of the first octant." << endl;
    cerr << "The connected component arising from a jump size of " << jumpSize << " is finite." << endl;
}