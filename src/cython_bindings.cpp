#include "../include/cython_bindings.hpp"
#include "../include/QuadrantSieve.hpp"
#include "../include/DonutSieve.hpp"
#include "../include/SegmentedSieve.hpp"
#include <cmath>


vector<pair<long, long>> gPrimes(long x) {
    // Show display if passed argument is large.
    bool display = x >= (long)pow(10, 7);
    DonutSieve s(x, display);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    vector<pair<long, long>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    return pairP;
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
