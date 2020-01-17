#include "../include/cython_bindings.hpp"
#include "../include/QuadrantSieve.hpp"
#include "OctantDonutSieve.hpp"
#include "../include/BlockSieve.hpp"
#include <iostream>
#include <cmath>


// Use this for ints?
vector<pair<uint32_t, uint32_t>> gPrimes(uint64_t x) {
    if (x >= pow(2, 30)) {
        cerr << "Norm bound too large for this method. Instead call gPrimesAsArray()" << endl;
        exit(1);
    }
    // Trigger verbose mode if passed argument is large.
    bool verbose = x >= (uint64_t)pow(10, 7);
    OctantDonutSieve s(x, verbose);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    vector<pair<uint32_t, uint32_t>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    cout << "Sending array over to python...." << endl;
    return pairP;
}


// And have this contain longs?
pair<uint32_t *, uint64_t> gPrimesAsArray(uint64_t x) {
    OctantDonutSieve s(x);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    // Creating a 1-dimensional array to hold big primes -- this way we can avoid array of pointers.
    uint64_t size = gintP.size();
    auto *P = new uint32_t[2 * size];  // declaring the array
    for (uint64_t i = 0; i < size; i++) {
        P[2 * i] = gintP[i].a;
        P[2 * i + 1] = gintP[i].b;
    }
    cout << P << endl;
    return pair<uint32_t *, uint64_t> {P, 2 * size};
}

vector<pair<uint32_t , uint32_t>> gPrimesSegment(uint32_t x, uint32_t y, uint32_t z) {
    // Show display if passed argument is large.
    bool display = z >= (uint32_t)pow(10, 5);
    BlockSieve s(x, y, z, display);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    vector<pair<uint32_t, uint32_t>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    return pairP;
}

uint64_t gPrimesCount(uint64_t x) {
    // Show display if passed argument is large.
    bool display = x >= (uint64_t)pow(10, 9);
    OctantDonutSieve s(x, display);
    s.run();
    return s.getCountBigPrimes();
}

uint64_t gPrimesSegmentCount(uint32_t x, uint32_t y, uint32_t z) {
    // Show display if passed argument is large.
    bool display = z >= (uint32_t)pow(10, 5);
    BlockSieve s(x, y, z, display);
    s.run();
    return s.getCountBigPrimes();
}
