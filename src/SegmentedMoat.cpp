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
uint32_t SegmentedMoat::previousdy;
uint64_t SegmentedMoat::blockSize;
uint64_t SegmentedMoat::sievingPrimesNormBound;
vector<gint> SegmentedMoat::sievingPrimes;
vector<gint> SegmentedMoat::nearestNeighbors;
vector<vector<gint>> SegmentedMoat::leftBoundary;
vector<uint64_t> SegmentedMoat::componentSizes;


// Call this static setter method before any instances of this class are created.
void SegmentedMoat::setStatics(double js, bool vb) {
    if (vb) {
        cerr << "Setting static variables..." << endl;
    }
    verbose = vb;

    // using tolerance with jumpSize
    double tolerance = pow(10, -3);
    jumpSize = js + tolerance;

    if (jumpSize < 3) {
        cerr << "Jump size is too small; instead call OctantMoat." << endl;
        exit(1);
    } else if (jumpSize >= 6) {
        cerr << "Jump size is too large for this implementation!" << endl;
        exit(1);
    }

    // Ideally should align to cache size, but also needs to be large enough to
    // hold more than solely the boundary data. Below are some parameters that
    // have worked in tests.
    if (jumpSize < 4) {
        blockSize = pow(10, 6);
    } else if (jumpSize < 4.1) {
        blockSize = pow(10, 7);
    } else if (jumpSize < 4.45) {
        blockSize = pow(10, 8);
    } else {
        blockSize = pow(10, 9);
    }

    // Bound on the norm of pre-computed primes. Arbitrary initial value.
    sievingPrimesNormBound = uint64_t(pow(10, 4));
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
    previousdy = floor(jumpSize);

    // Initializing leftBoundary to contain only 2 + i. Note that we need
    // jumpSize to be at least sqrt(2) to actually reach this prime from the origin.
    // Because of parity constraints in nearestNeighbors, we do not initialize with
    // 1 + i, but we still count it in componentSizes.
    vector<gint> component = {gint(2, 1)};
    // Initializing componentSizes
    componentSizes.push_back(2);  // primes 1 + i and 2 + i

    // Throwing in more gints and updating componentSizes if jumpSize is larger
    if (jumpSize - 1 > 3) {
        component.emplace_back(3, 0);
        component.emplace_back(3, 2);
        componentSizes[0] += 2;
    }
    if (jumpSize - 1 > 4) {
        component.emplace_back(4, 1);
        componentSizes[0]++;
    }
    if (jumpSize - 1 > 5) {
        component.emplace_back(5, 2);
        component.emplace_back(5, 4);
        componentSizes[0] += 2;
    }
    leftBoundary.push_back(component);
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
        , hasComponentPropagated(componentSizes.size(), false)  // no component propagated yet
{
    if (verbose) {
         cerr << "\n##################################################################" << endl;
         cerr << "Working within block having lower left corner at: " << x << " " << 0 << endl;
         cerr << "Size of block is: " << dx << " " << dy  << "\n" << endl;
    }

    if (2 * jumpSize > dx) {
        cerr << "Blocksize not large enough to fit both boundaries within sieveArray!" << endl;
        exit(1);
    }

    // Setting the 2D array holding component IDs within leftBoundary
    for (uint32_t a = 0; a < jumpSize - 1; a++) {
        vector<uint32_t> column(previousdy, 0);  // should be previous dy
        leftComponentLookUp.push_back(column);
    }
    for (uint32_t index = 0; index < leftBoundary.size(); index++) {
        for (gint g : leftBoundary[index]) {
            leftComponentLookUp[g.a][g.b] = index;
        }
        // Pushing empty vector onto rightBoundary for each distinct component
        vector<gint> component;
        rightBoundary.push_back(component);
    }
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


// Convention is to explore gint with smaller component index first so that
// merges go from large index into small index.
void SegmentedMoat::exploreComponent(uint32_t startingIndex, bool startingFromLeft) {
    uint64_t count = 0;  // everything in leftBoundary has been counted previously
    vector<gint> toExplore;
    if (startingFromLeft) {
        toExplore = leftBoundary[startingIndex];
    } else {  // removing singleton from rightBoundary since we'll push it back
        toExplore = rightBoundary[startingIndex];
        rightBoundary[startingIndex].clear();
    }

    // Indicating with sieveArray that these gints have been visited
    for (gint g : toExplore) { sieveArray[g.a][g.b] = false; }

    do {
        gint p = toExplore.back();
        toExplore.pop_back();

        // p is inside leftBoundary, so it was already counted in last iteration
        if (p.a < jumpSize - 1) {
            uint32_t index = leftComponentLookUp[p.a][p.b];
            if (index != startingIndex) {  // merging counts!
                // Adding the count from componentIndex to the current count,
                // then setting it to 0 to free it up for later use.
                count += componentSizes[index];
                componentSizes[index] = 0;

                // Pushing new component on toExplore, and updating
                // leftComponentLookUp. Because component indices are explored
                // in increasing order, none of the pushed gints below will have
                // been previously explored. We will end up exploring p twice;
                // p is only counted on the initial merge.
                for (gint g : leftBoundary[index]) {
                    leftComponentLookUp[g.a][g.b] = startingIndex;
                    toExplore.push_back(g);
                    sieveArray[g.a][g.b] = false;
                }
            }
        } else {  // haven't been to p in previous iteration
            count++;
            // p has punched through or started within the right boundary
            if (p.a >= floor(dx - jumpSize + 1)) {
                hasComponentPropagated[startingIndex] = true;
                rightBoundary[startingIndex].push_back(p);
            }
        }

        for (const gint &q : nearestNeighbors) {
            gint h = p + q;
            if (h.a >= 0 && h.a < dx && h.b >= 0 && h.b < dy && h.b <= x + h.a && sieveArray[h.a][h.b]) {
                toExplore.push_back(h);
                sieveArray[h.a][h.b] = false;  // indicating a visit here so we don't push back h again
            }
        }
    } while (!toExplore.empty());

    // Updating component count.
    componentSizes[startingIndex] += count;
}


void SegmentedMoat::exploreLeftBoundary() {
    // Stepping through the components in increasing order
    for (uint32_t index = 0; index < leftBoundary.size(); index++) {
        if (componentSizes[index]) {  // nonempty
            exploreComponent(index);
        }
    }

    // Now cleaning up components with index > 0 which have not propagated
    for (uint32_t index = 1; index < componentSizes.size(); index++) {
        if (!hasComponentPropagated[index]) {
            componentSizes[index] = 0;
        }
    }
}


void SegmentedMoat::exploreRightBoundary() {
    for (uint32_t a = floor(dx - jumpSize + 1); a < dx; a++) {
        for (uint32_t b = 0; b < dy; b++) {
            // Checking if unvisited, prime, and within first octant.
            if (sieveArray[a][b] && b < a + x) {
                // Finding an available index, or pushing new one.
                uint32_t index = 1;
                for (; index < componentSizes.size(); index++) {
                    if (componentSizes[index] == 0) {
                        // We found available index somewhere within componentCounts
                        // This will now be the index of this newly discovered component
                        break;
                    }
                }
                // Made it through componentCounts without finding an available
                // index. Creating a new one.
                if (index == componentSizes.size()) {
                    componentSizes.push_back(0);
                    hasComponentPropagated.push_back(true);
                    vector<gint> component;  // pushing empty component
                    rightBoundary.push_back(component);
                }
                rightBoundary[index].emplace_back(a, b);
                exploreComponent(index, false);
            }
        }
    }
}

void SegmentedMoat::runSegment() {
    exploreLeftBoundary();
    if (hasMainComponentPropagated()) {  // allows for early exit
        exploreRightBoundary();

        // Updating the static variable leftBoundary for next iteration
        leftBoundary = rightBoundary;
        for (auto &component : leftBoundary) {
            for (auto &g : component) {
                g.a -= floor(dx - jumpSize + 1);  // modifying gint
            }
        }
        // Updating static variable previousdy
        previousdy = dy;
    }
}


bool SegmentedMoat::hasMainComponentPropagated() {
    return hasComponentPropagated[0];
}


// Call this function after setStatics() has been run.
uint64_t SegmentedMoat::getCountMainComponent() {
    uint32_t x = 0;  // lower left corner of current block
    bool hasMainComponentPropagated;

    do {
        // Updating parameters dx and dy to pass to instance of SegmentedMoat.
        // Want: dx * dy = blockSize.
        // Also need: dy = x + dx so that next block goes all the way up to line y = x in complex plane.
        // Eliminating dy, these two equations give a quadratic in dx.
        uint32_t dx = floor(sqrt(blockSize + double(x) * double(x) / 4.0) - double(x) / 2.0);
        uint32_t dy = x + dx;

        // calling instance
        SegmentedMoat s(x, dx, dy);
        s.callSieve();
        s.runSegment();
        hasMainComponentPropagated = s.hasMainComponentPropagated();

        // updating x for next iteration
        x += floor(dx - jumpSize + 1);
    } while (hasMainComponentPropagated);
    return componentSizes[0];
}