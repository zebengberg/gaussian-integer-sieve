#include <iostream>
#include <random>
#include "../include/OctantSieve.hpp"
#include "../include/OctantDonutSieve.hpp"
#include "../include/BlockSieve.hpp"
#include "../include/BlockDonutSieve.hpp"
#include "../include/SectorSieve.hpp"
#include "../include/Moat.hpp"
using namespace std;


int main() {
    cout << "\nTesting and Timing Octant Sieve and Octant Donut Sieve\n" << endl;
    cout << " | norm bound | # of primes including associates | octant sieve time | octant donut sieve time | " << endl;
    cout << " |------------|----------------------------------|-------------------|-------------------------| " << endl;

    for (int j = 20; j <= 30; j++) {
        auto startTime = chrono::high_resolution_clock::now();
        OctantSieve o(pow(2, j), false);
        o.run();
        vector<gint> oP = o.getBigPrimes(false);  // not sorting yet
        auto endTime = chrono::high_resolution_clock::now();
        auto totalTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);
        double octantTime = double(totalTime.count()) / 1000.0;

        startTime = chrono::high_resolution_clock::now();
        OctantDonutSieve d(pow(2, j), false);
        d.run();
        vector<gint> dP = d.getBigPrimes(false);  // not sorting yet
        endTime = chrono::high_resolution_clock::now();
        totalTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);
        double donutTime = double(totalTime.count()) / 1000.0;

        cout << " | 2^" << j
             << " | " << 4 * oP.size()
             << " | " << octantTime
             << " s | " << donutTime
             << " s | " << endl;

        // Sorting generated primes and checking if two lists are equal.
        sort(oP.begin(), oP.end());
        sort(dP.begin(), dP.end());
        assert(oP == dP);
    }

    cout << "\nTesting and Timing Block Sieve and Block Donut Sieve\n" << endl;
    cout << " | block | # of primes | block sieve time | block donut sieve time | " << endl;
    cout << " |-------|-------------|------------------|------------------------| " << endl;

    random_device rd;
    uniform_int_distribution<int> distInt(1, 9);
    for (int j = 1; j <= 10; j++) {
        // Arbitrary testing parameters.
        uint32_t x = distInt(rd) * pow(10, 7);
        uint32_t y = distInt(rd) * pow(10, 7);
        uint32_t dx = distInt(rd) * 100;
        uint32_t dy = distInt(rd) * 100;

        auto startTime = chrono::high_resolution_clock::now();
        BlockSieve b(x, y, dx, dy, false);
        b.run();
        vector<gint> bP = b.getBigPrimes(false);  // not sorting yet
        auto endTime = chrono::high_resolution_clock::now();
        auto totalTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);
        double blockTime = double(totalTime.count()) / 1000.0;

        startTime = chrono::high_resolution_clock::now();
        BlockDonutSieve d(x, y, dx, dy, false);
        d.run();
        vector<gint> dP = d.getBigPrimes(false);  // not sorting yet
        endTime = chrono::high_resolution_clock::now();
        totalTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);
        double donutTime = double(totalTime.count()) / 1000.0;

        cout << " | [" << x << ", " << x + dx
             << ") x [" << y << ", " << y + dy
             << ") | " << bP.size()
             << " | " << blockTime
             << " s | " << donutTime
             << " s | " << endl;

        // Sorting generated primes and checking if two lists are equal.
        sort(bP.begin(), bP.end());
        sort(dP.begin(), dP.end());
        assert(bP == dP);
    }

    cout << "\nTesting and Timing Sector Sieve\n" << endl;
    cout << " | alpha | beta | beta - alpha | norm bound | # of primes | time | " << endl;
    cout << " |-------|------|--------------|------------|-------------|------| " << endl;
    uniform_real_distribution<long double> distReal(0.0, M_PI_4);
    for (int j = 20; j <= 30; j++) {
        uint64_t x = pow(2, j);
        long double alpha = distReal(rd);
        long double beta = distReal(rd);
        if (beta < alpha) {
            long double temp = beta;
            beta = alpha;
            alpha = temp;
        }

        auto startTime = chrono::high_resolution_clock::now();
        SectorSieve s(x, alpha, beta, false);
        s.run();
        vector<gint> sP = s.getBigPrimes(false);  // not sorting yet
        auto endTime = chrono::high_resolution_clock::now();
        auto totalTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);
        double sectorTime = double(totalTime.count()) / 1000.0;
        cout << " | " << alpha
             << " | " << beta
             << " | " << beta - alpha
             << " | 2^" << j
             << " | " << sP.size()
             << " | " << sectorTime
             << " | " << endl;

        OctantDonutSieve d(x, false);
        d.run();
        vector<gint> dP = d.getBigPrimes();
        // Removing gints outside of sector
        auto res = remove_if(dP.begin(), dP.end(), [&](gint g) { return (g.arg() >= beta) || (g.arg() < alpha); });
        dP.erase(res, dP.end());

        // Sorting generated primes and checking if two lists are equal.
        sort(sP.begin(), sP.end());
        sort(dP.begin(), dP.end());
        assert(sP == dP);
    }

    cout << "\nTiming Thin Sectors\n" << endl;
    cout << " | alpha | beta | beta - alpha | norm bound | # of primes | time | " << endl;
    cout << " |-------|------|--------------|------------|-------------|------| " << endl;
    for (int j = 30; j <= 40; j++) {
        uint64_t x = pow(2, j);
        double alpha = distReal(rd);
        double beta = alpha + pow(2, 20 - j);

        auto startTime = chrono::high_resolution_clock::now();
        SectorSieve s(x, alpha, beta, false);
        s.run();
        vector<gint> sP = s.getBigPrimes(false);  // not sorting yet
        auto endTime = chrono::high_resolution_clock::now();
        auto totalTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);
        double sectorTime = double(totalTime.count()) / 1000.0;
        cout << " | " << alpha
             << " | " << beta
             << " | 2^" << 20 - j
             << " | 2^" << j
             << " | " << sP.size()
             << " | " << sectorTime
             << " | " << endl;
    }

    cout << "Exploring the Gaussian moat." << endl;
    // Data taken from the 1998 paper "A stroll through the Gaussian primes".

    OctantMoat m(10000, 2);
    m.exploreComponent(0, 0);
    assert(m.getComponentSize() == 92);
    assert(m.getComponentMaxElement() == gint(42, 17));

    m = OctantMoat(100000, 3);
    m.exploreComponent(0, 0);
    assert(m.getComponentSize() == 380);
    assert(m.getComponentMaxElement() == gint(84, 41));

    m = OctantMoat(1000000, 3.2);
    m.exploreComponent(0, 0);
    assert(m.getComponentSize() == 31221);
    assert(m.getComponentMaxElement() == gint(976, 311));

    m = OctantMoat(10000000, 4);
    m.exploreComponent(0, 0);
    assert(m.getComponentSize() == 347638);
    assert(m.getComponentMaxElement() == gint(3297, 2780));

    m = OctantMoat(10000000, 4.3);
    m.exploreComponent(0, 0);
    assert(m.getComponentSize() == 2386129);
    assert(m.getComponentMaxElement() == gint(8174, 6981));

    return 0;
}
