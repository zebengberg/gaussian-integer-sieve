// A library of functions to be wrapped into a python module.

#include "../include/cython_bindings.hpp"
#include "../include/OctantDonutSieve.hpp"
#include "../include/SectorSieve.hpp"
#include "../include/Moat.hpp"
#include <iostream>
#include <cmath>
#include <numeric>


// Creating a 1-dimensional arrays to hold big primes; this way we can avoid
// an array of pointers which might be needed for 2d array. This will be fed
// into a memory view object a la cython to give numpy access. The keyword
// "new" will prevent the array from decaying into garbage.
// This function is called in various functions below.
// In the cython file gintsieve.pyx, there is an inverse function to this one.
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

// Counting Gaussian primes and associates upto a given norm.
uint64_t gPrimesToNormCount(uint64_t x) {
    if (x < 2) {
        return 0;
    } else if (x < 5) {
        return 4;
    } else {
        // Show display if passed argument is large.
        bool verbose = x >= (uint64_t) pow(10, 9);
        OctantDonutSieve s(x, verbose);
        s.run();
        return s.getCountBigPrimes();
    }
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
    vector<gint> gintP;
    if (x >= 5) {
        // Show display if passed argument is large.
        bool verbose = x >= (uint64_t)pow(10, 9);
        OctantDonutSieve s(x, verbose);
        s.run();
        gintP = s.getBigPrimes();
    } else if (x >= 2) {
        gintP.emplace_back(1, 1);
    }
    return gintVectorToArray(gintP);
}

// Getting Gaussian primes in sector upto given norm. Return pointer to an array.
pair<int32_t *, uint64_t> gPrimesInSectorAsArray(uint64_t x, double alpha, double beta) {
    // Show display if passed argument is large.
    bool verbose = x >= (uint64_t)pow(10, 9);
    SectorSieve s(x, alpha, beta, verbose);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    return gintVectorToArray(gintP);
}

// Getting Gaussian primes in block. Return pointer to an array.
pair<int32_t *, uint64_t> gPrimesInBlockAsArray(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy) {
    // Show display if passed argument is large.
    bool verbose = dx * dy >= (uint64_t)pow(10, 9);
    BlockSieve s(x, y, dx, dy, verbose);
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


// Wrapper functions to access moat data
pair<int32_t *, uint64_t> moatMainComponent(double jumpSize) {
    OctantMoat m(jumpSize);
    m.exploreComponent(0, 0);
    vector<gint> component = m.getCurrentComponent();
    // pushing gint with max norm onto end of vector
    gint g = m.getComponentMaxElement();
    component.push_back(g);
    return gintVectorToArray(component);
}

vector<pair<int32_t *, uint64_t>> moatComponentsToNorm(double jumpSize, uint64_t x) {
    OctantMoat m(jumpSize, x);
    m.exploreAllComponents();
    vector<vector<gint>> allComponents = m.getAllComponents();
    vector<pair<int32_t *, uint64_t>> toReturn;
    toReturn.reserve(allComponents.size());  // pre-allocating size
    for (const vector<gint>& v : allComponents) {
        toReturn.push_back(gintVectorToArray(v));
    }
    return toReturn;
}

vector<pair<int32_t *, uint64_t>> moatComponentsInBlock(double jumpSize,
        int32_t x, int32_t y, int32_t dx, int32_t dy) {
    BlockMoat m(jumpSize, x, y, dx, dy);
    m.run();  // from parent BlockSieve
    m.exploreAllComponents();

    vector<vector<gint>> allComponents = m.getAllComponents();
    vector<gint> edges = m.getEdges();
    // putting edges in allComponents while we pass data to python
    allComponents.push_back(edges);

    vector<pair<int32_t *, uint64_t>> toReturn;
    toReturn.reserve(allComponents.size());  // pre-allocating size
    for (const vector<gint>& v : allComponents) {
        toReturn.push_back(gintVectorToArray(v));
    }
    return toReturn;
}


// Constructor
BlockMoat::BlockMoat(double js, uint32_t x, uint32_t y, uint32_t dx, uint32_t dy)
// Calling BlockSieve's constructor
        : BlockSieve(x, y, dx, dy, false) // not letting this be verbose
        , jumpSize(js)
        , dx(dx)
        , dy(dy)
{
    jumpSize += pow(10, -3);  // adding a small tolerance
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
}

void BlockMoat::exploreComponent(gint g0) {
    // reset current component
    currentComponent.clear();
    currentComponent.push_back(g0);

    // Because of parity constraint in nearestNeighbors, this method will
    // never visit the ramifying prime 1 + i.
    vector<gint> toExplore;
    toExplore.push_back(g0);

    do {
        gint p = toExplore.back();
        toExplore.pop_back();
        sieveArray[p.a][p.b] = false;  // indicating a visit
        for (const gint &q : nearestNeighbors) {
            gint g = p + q;

            // Pushing neighbors onto the vector toExplore
            if (g.a >= 0 && g.a < dx && g.b >= 0 && g.b < dy && sieveArray[g.a][g.b]) {
                currentComponent.push_back(g);
                toExplore.push_back(g);
                edges.push_back(p);
                edges.push_back(g);
            }
        }
    }
    while (!toExplore.empty());
}

void BlockMoat::exploreAllComponents() {
    for (int32_t a = 0; a < dx; a++) {
        for (int32_t b = 0; b < dy; b++) {
            if (sieveArray[a][b]) {
                gint g(a, b);
                exploreComponent(g);
                allComponents.push_back(currentComponent);
            }
        }
    }
}

vector<vector<gint>> BlockMoat::getAllComponents() { return allComponents; }

vector<gint> BlockMoat::getEdges() { return edges; }