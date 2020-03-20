// Function to count size of connected component at origin using many blocks
// rather than a single gargantuan call to OctantMoat.

// ALGORITHM:
// Break the first octant into many trapezoidal vertical strips. Within each
// vertical strip, instantiate a BlockMoat-type object to explore all components
// that eventually meander to the right wall of the block. Along this right
// boundary, keep track of the size of the component as well as an ID for the
// component ending there. Because the boundary is closer to a 1-dimensional
// geometric object rather than the 2-dimensional block, tracking component
// status along this boundary will require less storage. Specifically, this
// algorithm will only require O(sqrt(X)) storage to find the size of the
// component containing the origin out to points with norm X.
//
// Exploring components in the next block, we simply start at the previous
// boundary. We only care about persistent components that propagate to the
// right. Islands and lakes will be forgotten with this algorithm; we only seek
// to find the size of the component containing the origin. Counts will be
// merged when components come together in future blocks.

/*
 * Three parts to algorithm:
 * 1. For each gint g in left boundary with component 1, explore at g. When encounter
 * other gints in leftBoundary with different components, merge counts. When cross into
 * rightBoundary, mark off visits.
 * 2. Go through each gint g in left boundary with component != 1. Follow procedure as above.
 * 3. For each unvisited prime g in the right boundary, explore there to update the count.
 *
 */

#include <iostream>
#include "../include/Moat.hpp"
#include "../include/OctantDonutSieve.hpp"

// Need to first declare static member variables here in the cpp source.
bool SegmentedMoat::verbose;
double SegmentedMoat::jumpSize;
uint64_t SegmentedMoat::sievingPrimesNormBound;
vector<gint> SegmentedMoat::sievingPrimes;
vector<gint> SegmentedMoat::nearestNeighbors;
vector<vector<uint32_t>> SegmentedMoat::leftBoundary;
vector<uint64_t> SegmentedMoat::componentCounts;
bool SegmentedMoat::mainComponentPunchedThrough;


// Call this static setter method before any instances of this class are created.
void SegmentedMoat::setStatics(double js, bool vb) {
    if (vb) {
        cerr << "Setting static variables..." << endl;
    }
    verbose = vb;

    jumpSize = js;
    // using tolerance with jumpSize
    double tolerance = pow(10, -3);
    jumpSize += tolerance;

    if (jumpSize < 3) {
        cerr << "Jump size is too small; instead call OctantMoat." << endl;
        exit(1);
    }

    // Bound on the norm of pre-computed primes. Initial value is arbitrary.
    sievingPrimesNormBound = uint64_t(pow(10, 8));
    setSievingPrimes();

    // Setting nearest neighbors.
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

    // Setting leftBoundary to only contain origin.
    for (uint32_t i = 0; i < jumpSize; i++) {
        vector<uint32_t> column(1, 0);  // ID of 0 indicates unvisited
        leftBoundary.push_back(column);
    }
    leftBoundary[0][0] = 1;  // ID of component at origin is 1

    // Initializing component counts vector.
    componentCounts.push_back(0);
    componentCounts.push_back(0);
}

void SegmentedMoat::setSievingPrimes() {
    if (verbose) {
        cerr << "Precomputing sieving primes." << endl;
    }
    OctantDonutSieve d(sievingPrimesNormBound);
    d.run();
    sievingPrimes = d.getBigPrimes();
}

SegmentedMoat::SegmentedMoat(int32_t x, int32_t dx, int32_t dy)
// Calling BlockSieve's constructor
        : BlockSieve(x, 0, dx, dy, false) // not letting this be verbose
        , x(x)
        , dx(dx)
        , dy(dy)
{
    if (verbose) {
        cerr << "Working within block having lower left corner at: " << x << " " << 0 << endl;
    }

    // Setting rightBoundary to contain jumpSize-number of columns; should be
    // indexed by last jumpSize columns of the sieveArray.
    for (uint32_t i = 0; i < jumpSize; i++) {
        vector<uint32_t> column(dy, 0);  // 0 indicates unvisited
        rightBoundary.push_back(column);
    }

    mainComponentPunchedThrough = false;
}


// Cannot call virtual methods of BlockSieve parent from BlockMoat constructor.
void SegmentedMoat::callSieve() {
    // Checking to make sure there are enough primes within sievingPrimes
    gint last_g = sievingPrimes.back();
    while (last_g.norm() < isqrt(maxNorm)) {
        cerr << "Not enough pre-computed primes in static variable sievingPrimes." << endl;
        cerr << "Doubling the norm bound on pre-computed primes and computing more..." << endl;
        sievingPrimesNormBound *= 2;
        setSievingPrimes();
        last_g = sievingPrimes.back();
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

// Here g should be a gint in the left boundary of the sieve array. Convention
// is to explore gint with smaller component index first so that merges go from
// large index into small index.
void SegmentedMoat::exploreAtGint(gint g, uint32_t startingComponentIndex) {
    // to the right boundary. If not, entire component count can be forgotten.
    uint64_t count = 0;

    vector<gint> toExplore;
    toExplore.push_back(g);

    do {
        gint p = toExplore.back();
        toExplore.pop_back();
        sieveArray[p.a][p.b] = false;  // indicating a visit
        if (p.a < jumpSize) {  // p is inside leftBoundary
            uint64_t componentIndex = leftBoundary[p.a][p.b];
            if (componentIndex != startingComponentIndex) {  // merge counts!
                count += componentCounts[componentIndex];
                // Setting count to 0 to avoid over count when we encounter other
                // gints with index componentIndex.
                componentCounts[componentIndex] = 0;
            }  // else, p is in the same component as g, and we've already counted contribution
        } else {  // haven't yet been to p
            count++;
            if (p.a >= floor(dy - jumpSize)) {  // p has punched through into right boundary!
                if (startingComponentIndex == 1) {  // main component has successfully propagated
                    mainComponentPunchedThrough = true;
                }
                rightBoundary[p.a - floor(dy - jumpSize)][p.b] = startingComponentIndex;
                count++;
            }
        }

        for (const gint &q : nearestNeighbors) {
            gint h = p + q;
            if (h.a >= 0 && h.a <= dx && h.b >= 0 && h.b <= dy && h.b <= h.a && sieveArray[h.a][h.b]) {
                toExplore.push_back(h);
            }
        }
    } while (!toExplore.empty());

    // Updating component count.
    componentCounts[startingComponentIndex] += count;
}

void SegmentedMoat::exploreLeftBoundary() {
    // Going from gints with low component index to high component index.
    // Need to step through leftBoundary several times; consider re-engineering
    // if this is too slow.
    for (uint32_t index = 0; index < componentCounts.size(); index++) {
        for (uint32_t a = 0; a < jumpSize; a++) {
            for (uint32_t b = 0; b < dy; b++) {
                if (leftBoundary[a][b] == index) {
                    exploreAtGint(gint(a, b), index);
                }
            }
        }
    }
}

void SegmentedMoat::exploreRightBoundary() {
    for (uint32_t a = 0; a < jumpSize; a++) {
        for (uint32_t b = 0; b < dy; b++) {
            uint32_t aAbsolute = floor(dy - jumpSize) + a;
            // Checking if unvisited, prime, and within first octant.
            if (rightBoundary[a][b] == 0 && sieveArray[aAbsolute][b] && b < aAbsolute) {
                // Finding an available opening in segmentedCounts, or pushing new one.
                uint32_t index = 1;
                for (; index < componentCounts.size(); index++) {
                    if (componentCounts[index] == 0) {
                        // We found available index somewhere within componentCounts
                        componentCounts[index]++;
                        break;
                    }
                }
                if (index == componentCounts.size()) {  // appending new entry to vector
                    componentCounts.push_back(1);
                }
                exploreAtGint(gint(aAbsolute, b), index);
            }
        }
    }
}

void SegmentedMoat::runSegment() {
    if (verbose) {
        cerr << "Exploring left boundary..." << endl;
    }
    exploreLeftBoundary();
    if (verbose) {
        cerr << "Exploring right boundary..." << endl;
    }
    exploreRightBoundary();
    leftBoundary = rightBoundary;
}

// Call this static method after setStatics() has been run.
uint64_t SegmentedMoat::countComponent() {
    uint32_t x = 0;  // lower left corner of current block
    // size of each block; ideally should align to cache size, but this also
    // needs to be large enough to hold more than boundary data
    // TODO: should blocksize be static?
    uint64_t blockSize = pow(10, 8);
    int32_t dx = pow(10, 4);
    int32_t dy = pow(10, 4);

    do {
        // calling instance
        SegmentedMoat s(x, dx, dy);
        s.callSieve();
        s.runSegment();

        // Updating parameters for next call.
        // Want: dx * dy = blockSize.
        // Also need: dy = x + dx so that next block goes all the way up to line y = x in complex plane.
        // Eliminating dy, these two equations give a quadratic in dx.
        x += dx - floor(jumpSize);
        dx = floor(sqrt(blockSize + x * x / 4.0) - x / 2.0);
        dy = x + dx;

    } while (mainComponentPunchedThrough);
    return componentCounts[1];
}