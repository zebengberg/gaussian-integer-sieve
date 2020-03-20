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





// Need to first declare static member variables here in source.
bool BlockMoat::verbose;
double BlockMoat::jumpSize;
int32_t BlockMoat::dx;
int32_t BlockMoat::dy;
uint64_t BlockMoat::sievingPrimesNormBound;
vector<gint> BlockMoat::sievingPrimes;
vector<gint> BlockMoat::nearestNeighbors;

// Call this static setter method before any instances of this class are created.
void BlockMoat::setStatics(int32_t realPart, double js, bool vb) {
    if (vb) {
        cerr << "Setting static variables..." << endl;
    }
    verbose = vb;
    jumpSize = js;
    dx = 1000;
    dy = 10000;
    // Bound on the norm of pre-computed primes. The factor 1.2 gives some
    // wiggle room in case there are many moves to the right.
    sievingPrimesNormBound = uint64_t(1.2 * (sqrt(2) * realPart + dx * dy));

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

    if (verbose) {
        cerr << "Precomputing sieving primes." << endl;
    }
    OctantDonutSieve d(sievingPrimesNormBound);
    d.run();
    sievingPrimes = d.getBigPrimes();
}

// Constructor
BlockMoat::BlockMoat(int32_t x, int32_t y)
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
    }
}

// Cannot call virtual methods of BlockSieve parent from BlockMoat constructor.
void BlockMoat::callSieve() {
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
bool BlockMoat::exploreAtGint(int32_t a, int32_t b, bool upperWallFlag) {
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
            if ((g.a >= 0) && (g.a < dx) && (g.b >= 0) && (g.b < dy) && sieveArray[g.a][g.b]) {
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
        if (verbose) {
            cerr << "Largest real-part reached starting from left hand wall: " << farthestRight << endl;
            cerr << "Number of visited primes: " << countVisited << "\n" << endl;
        }
        return {x, y + upperWallYPunch};
    }
}

// Function to use BlockMoat.
void verticalMoat(int32_t realPart, double jumpSize, bool verbose) {
    BlockMoat::setStatics(realPart, jumpSize, verbose);
    int32_t x = realPart;
    int32_t y = 0;
    int32_t consecutiveStepsRight = 0;

    // TODO: Update values of dx and dy depending on the size of the real parts reached
    while (y < x) {
        // TODO: Check if this object needs to be explicitly destructed
        BlockMoat b(x, y);
        b.callSieve();
        pair<int32_t, int32_t> p = b.getNextBlock();
        if (p.first != x) {
            consecutiveStepsRight++;
        } else {
            consecutiveStepsRight = 0;
        }
        if (consecutiveStepsRight > 10) {
            cerr << "Stepped right ten times in a row!" << endl;
            exit(1);
        }
        x = p.first;
        y = p.second;
    }
    cerr << "Gaussian moat present from real-axis to boundary of the first octant." << endl;
    cerr << "The connected component arising from a jump size of " << jumpSize << " is finite." << endl;
}




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

// Need to first declare static member variables here in the cpp source.
bool SegmentedMoat::verbose;
double SegmentedMoat::jumpSize;
uint64_t SegmentedMoat::sievingPrimesNormBound;
vector<gint> SegmentedMoat::sievingPrimes;
vector<gint> SegmentedMoat::nearestNeighbors;
vector<vector<uint64_t>> SegmentedMoat::leftBoundary;
vector<uint64_t> SegmentedMoat::componentCounts;

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
        vector<uint64_t> column(1, -1);  // -1 indicates unvisited
        leftBoundary.push_back(column);
    }
    leftBoundary[0][0] = 1;  // ID of component at origin is 1
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
        cerr << "Working within trapezoid having lower left corner at: " << x << " " << 0 << endl;
    }

    // Setting rightBoundary to contain jumpSize-number of columns; should be
    // indexed by last jumpSize columns of the sieveArray.
    for (uint32_t i = floor(dy - jumpSize); i < dx; i++) {
        vector<uint64_t> column(dy, 0);  // 0 indicates unvisited
        rightBoundary.push_back(column);
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

// Here g should be a gint in the left boundary of the sieve array. Convention
// is to explore gint with smaller component index first so that merges go from
// large index into small index.
void SegmentedMoat::exploreAtGint(gint g, uint64_t startingComponentIndex) {
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
            if (p.a >= floor(dy - jumpSize)) {  // p is actually inside rightBoundary
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
};

void SegmentedMoat::exploreLeftBoundary() {
    // Going from gints with low component index to high component index.
    // Need to step through leftBoundary several times; consider re-engineering
    // if this is too slow.
    for (uint64_t index = 0; index < componentCounts.size(); index++) {
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
    // TODO: write this
}


/*
 * Three parts to algorithm:
 * 1. For each gint g in left boundary with component 1, explore at g. When encounter
 * other gints in leftBoundary with different components, merge counts. When cross into
 * rightBoundary, mark off visits.
 * 2. Go through each gint g in left boundary with component != 1. Follow procedure as above.
 * 3. For each unvisited prime g in the right boundary, explore there to update the count.
 *
 */


// TODO: put different classes in different files; keep same header.













// Something old -- possibly delete.

uint64_t segmentedMoat(double jumpSize, bool verbose) {
    uint64_t count = 0;
    uint32_t realPart = 0;  // lower left corner of current block
    // size of each block; ideally should align to cache size, but this also
    // needs to be large enough to hold more than boundary data
    uint64_t blockSize = pow(10, 8);
    int32_t dx = pow(10, 4);
    int32_t dy = pow(10, 4);




    BlockMoat b(realPart, 0);
    b.callSieve();

    vector<pair<vector<gint>, uint64_t>> oldBoundary;  // list of boundary gints, count
    vector<pair<vector<gint>, uint64_t>> newBoundary;
    for (const auto& component : oldBoundary) {
        for (gint g : component.first) {
            b.exploreAtGint(g.a)
        }
    }



}