#include "../include/cython_bindings.hpp"
#include "../include/QuadrantSieve.hpp"
#include "../include/DonutSieve.hpp"
#include "../include/SegmentedSieve.hpp"
#include <iostream>
#include <cmath>


// Use this for ints?
vector<pair<long, long>> gPrimes(long x) {
    if (x >= pow(2, 32)) { // between 10^9 and 10^10
        cerr << "Norm bound too large for this method. Instead call gPrimes2dArray()" << endl;
        exit(1);
    }
    // Trigger verbose mode if passed argument is large.
    bool verbose = x >= (long)pow(10, 7);
    DonutSieve s(x, verbose);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    vector<pair<long, long>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    cout << "Sending array over to python...." << endl;
    return pairP;
}


// And have this contain longs?
pair<long *, unsigned long> gPrimesArray(long x) {
    DonutSieve s(x);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    // Creating a 1-dimensional array to hold big primes -- this way we can avoid array of pointers.
    unsigned long size = gintP.size();
    long *P = new long[2 * size];
    for (unsigned long i = 0; i < size; i++) {
        P[2 * i] = gintP[i].a;
        P[2 * i + 1] = gintP[i].b;
    }
    cout << P << endl;
    return pair<long *, unsigned long> {P, 2 * size};
}

vector<pair<long, long>> gPrimesSegment(long x, long y, long z) {
    // Show display if passed argument is large.
    bool display = z >= (long)pow(10, 5);
    SegmentedSieve s(x, y, z, display);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    vector<pair<long, long>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    return pairP;
}

unsigned long gPrimesCount(long x) {
    // Show display if passed argument is large.
    bool display = x >= (long)pow(10, 9);
    DonutSieve s(x, display);
    s.run();
    return s.getCountBigPrimes();
}

unsigned long gPrimesSegmentCount(long x, long y, long z) {
    // Show display if passed argument is large.
    bool display = z >= (long)pow(10, 5);
    SegmentedSieve s(x, y, z, display);
    s.run();
    return s.getCountBigPrimes();
}
