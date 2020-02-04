// A library of functions to be wrapped into a python module.

#include "../include/cython_bindings.hpp"
#include "../include/OctantDonutSieve.hpp"
#include "../include/SectorSieve.hpp"
#include "../include/BlockSieve.hpp"
#include "../include/OctantSieve.hpp"
#include <iostream>
#include <cmath>
#include <numeric>


// Creating a 1-dimensional arrays to hold big primes; this way we can avoid
// an array of pointers which might be needed for 2d array. This will be fed
// into a memory view object a la cython to give numpy access. The keyword
// "new" will prevent the array from decaying into garbage.
// This function is called in various functions below.
pair<int32_t *, uint64_t> gintVectorToArray(vector<gint> v) {
    // Creating a 1-dimensional array to hold big primes; this way we can avoid
    // an array of pointers which might be needed for 2d array.
    uint64_t size = v.size();
    auto *a = new int32_t[2 * size];  // declaring the array
    for (uint64_t i = 0; i < size; i++) {
        a[2 * i] = v[i].a;
        a[2 * i + 1] = v[i].b;
    }
    return pair<int32_t *, uint64_t> {a, 2 * size};
}


// Getting Gaussian primes upto a given norm. Python will hold these as a list
// of tuples. Each tuples holds the real and imaginary parts of each prime.
vector<pair<int32_t, int32_t>> gPrimesToNorm(uint64_t x) {
    if (x >= pow(2, 30)) {
        cerr << "Large norm bound makes copying C++ containers to python\n"
             << "containers too slow. Call gPrimesAsArray() instead.\n" << endl;
        exit(1);
    }
    // Trigger verbose mode if passed argument is large.
    bool verbose = x >= (uint64_t)pow(10, 7);
    OctantDonutSieve s(x, verbose);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    vector<pair<int32_t, int32_t>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    cerr << "Copying sieve array to python list of tuples...." << endl;
    return pairP;
}

//  Getting Gaussian primes in sector upto a given norm.
vector<pair<int32_t, int32_t>> gPrimesInSector(uint64_t x, double alpha, double beta) {
    if (x >= pow(2, 30)) {
        cerr << "Large norm bound makes copying C++ containers to python\n"
             << "containers too slow. Call gPrimesAsArray() instead.\n" << endl;
        exit(1);
    }
    // Trigger verbose mode if passed argument is large.
    bool verbose = x >= (uint64_t)pow(10, 7);
    SectorSieve s(x, alpha, beta, verbose);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    vector<pair<int32_t, int32_t>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    cerr << "Copying sieve array to python list of tuples...." << endl;
    return pairP;
}

// Gaussian integers in block defined by intervals [x, x + dx) and [y, y + dy).
vector<pair<int32_t, int32_t>> gPrimesInBlock(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy) {
    if (dx * dy >= pow(2, 30)) {
        cerr << "Large norm bound makes copying C++ containers to python\n"
             << "containers too slow. Call gPrimesAsArray() instead.\n" << endl;
        exit(1);
    }
    // Trigger verbose mode if passed argument is large.
    bool verbose = dx * dy >= (uint64_t)pow(10, 7);
    BlockSieve s(x, y, dx, dy, verbose);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    vector<pair<int32_t, int32_t>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    cerr << "Copying sieve array to python list of tuples...." << endl;
    return pairP;
}

// Counting Gaussian primes and associates upto a given norm.
uint64_t gPrimesToNormCount(uint64_t x) {
    // Show display if passed argument is large.
    bool verbose = x >= (uint64_t)pow(10, 9);
    OctantDonutSieve s(x, verbose);
    s.run();
    return s.getCountBigPrimes();
}

// Counting Gaussian primes in sector upto a given norm.
uint64_t gPrimesInSectorCount(uint64_t x, double alpha, double beta) {
    // Show display if passed argument is large.
    bool verbose = x >= (uint64_t)pow(10, 9);
    SectorSieve s(x, alpha, beta, verbose);
    s.run();
    return s.getCountBigPrimes();
}

// Counting Gaussian primes in sector upto a given norm.
uint64_t gPrimesInBlockCount(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy) {
    // Show display if passed argument is large.
    bool verbose = dx * dy >= (uint64_t)pow(10, 9);
    BlockSieve s(x, y, dx, dy, verbose);
    s.run();
    return s.getCountBigPrimes();
}

// Getting Gaussian primes upto a given norm. Passing a pointer to an array in
// c++, so data is not explicitly copied between memory controlled by c++ and
// python.
pair<int32_t *, uint64_t> gPrimesToNormAsArray(uint64_t x) {
    OctantDonutSieve s(x);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    return gintVectorToArray(gintP);
}

// Getting Gaussian primes in sector upto given norm. Return pointer to an array.
pair<int32_t *, uint64_t> gPrimesInSectorAsArray(uint64_t x, double alpha, double beta) {
    SectorSieve s(x, alpha, beta);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    return gintVectorToArray(gintP);
}

// Getting Gaussian primes in block. Return pointer to an array.
pair<int32_t *, uint64_t> gPrimesInBlockAsArray(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy) {
    BlockSieve s(x, y, dx, dy);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    return gintVectorToArray(gintP);
}

// Getting statistics on the angular distribution of Gaussian primes to norm.
vector<uint64_t> angularDistribution(uint64_t x, uint32_t nSectors) {
    // A vector holding nSectors, each with value 0.
    vector<uint64_t> sectors(nSectors, 0);
    OctantDonutSieve s(x);
    s.run();
    vector<gint> gintP = s.getBigPrimes(false);  // not sorting
    cerr << "Putting primes into bins according to their angle...." << endl;
    for (gint g : gintP) {
        double angle = g.arg();
        auto sector = uint32_t(nSectors * angle / M_PI_4);
        // Only considering first octant (not second).
        // Inert primes cause outliers in first sector count;
        // perhaps should count with multiplicity 1/2.
        if (sector < nSectors) {
            sectors[sector]++;
        }
    }
    return sectors;
}

// Public methods in SectorRace class.
SectorRace::SectorRace(uint64_t x, uint64_t nBins, long double alpha, long double beta, long double gamma, long double delta)
    // initializer list
    : x(x)
    , nBins(nBins)
{
    cerr << "Running Sector Sieves...\n" << endl;
    SectorSieve s1(x, alpha, beta);
    SectorSieve s2(x, gamma, delta);
    s1.run();
    s2.run();
    // Not sorting big primes.
    firstSector = s1.getBigPrimes(false);
    secondSector = s2.getBigPrimes(false);
    setNormData();
}

void SectorRace::setNormData() {
    cerr << "Collecting norm data for sector race..." << endl;
    // Collecting norm data by creating an vector of all 0s with length nBins
    for (uint64_t i = 0; i < nBins; i++) {
        normData.push_back(0);
    }
    // Incrementing for gints in first sector; decrementing for gints in second sector.
    for (gint g : firstSector) {
        normData[g.norm() * nBins / x]++;
    }
    for (gint g : secondSector) {
        normData[g.norm() * nBins / x]--;
    }
    // Taking cumulative sum and overwriting normData with it.
    partial_sum(normData.begin(), normData.end(), normData.begin());
}

pair<int32_t *, uint64_t> SectorRace::getFirstSector() {
    return gintVectorToArray(firstSector);
}

pair<int32_t *, uint64_t> SectorRace::getSecondSector() {
    return gintVectorToArray(secondSector);
}

// Putting normData contents into an array; see comment above on memory view.
int32_t * SectorRace::getNormData() {
    auto *data = new int32_t[nBins];  // declaring the array
    for (uint32_t i = 0; i < nBins; i++) {
        data[i] = normData[i];
    }
    return data;
}


// Public methods in OctantMoat class.
OctantMoat::OctantMoat(uint64_t x, double jumpSize) : x(x), jumpSize(jumpSize)
{
    setSieveArray();
    setNeighbors();
    setToExplore();
    explore();
}

void OctantMoat::setSieveArray() {
    OctantSieve o(x);
    o.run();
    sieveArray = o.getSieveArray();
}

void OctantMoat::setNeighbors() {
    for (int32_t u = -int32_t(jumpSize); u < jumpSize; u++) {
        for (int32_t v = -int32_t(jumpSize); v < jumpSize; v++) {
            if (u * u + v * v <= jumpSize * jumpSize) {
                neighbors.emplace_back(gint(u, v));
            }
        }
    }
}

void OctantMoat::setToExplore() {
    toExplore.emplace_back(gint(0, 0));
}

void OctantMoat::explore() {
    uint32_t count = 0;
    while (!toExplore.empty()) {
        gint p = toExplore.back();
        toExplore.pop_back();
        for (gint q : neighbors) {
            gint g = p + q;
            // Checking if inside first octant and prime
            if ((g.a >= 0) && (g.b >= 0) && (g.b <= g.a) && sieveArray[g.a][g.b]) {
                toExplore.push_back(g);
                sieveArray[g.a][g.b] = false;
            }
        }
        explored.push_back(p);
        count++;
        if (count % 10 == 0) {
            cerr << '.';
        }
        if (count % 1000 == 0) {
            cerr << endl;
        }
    }
}

pair<int32_t *, uint64_t> OctantMoat::getExplored() {
    return gintVectorToArray(explored);
}

pair<int32_t *, uint64_t> OctantMoat::getUnexplored() {
    vector<gint> unexplored;
    for (uint32_t u = 0; u < sieveArray.size(); u++) {
        for (uint32_t v = 0; v < sieveArray[u].size(); v++) {
            if (sieveArray[u][v]) {
                unexplored.emplace_back(u, v);
            }
        }
    }
    return gintVectorToArray(unexplored);
}
