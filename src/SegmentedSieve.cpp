#include <iostream>
#include <cmath>
#include "../include/SegmentedSieve.hpp"
#include "../include/OctantSieve.hpp"
using namespace std;

// Using an initializer list
SegmentedSieve::SegmentedSieve(uint32_t x, uint32_t y, uint32_t z, uint32_t w, bool verbose)
    // Calling SieveTemplate constructor to set maxNorm
    : SieveTemplate<bool>((uint64_t)(x + w - 1) * (uint64_t)(y + z - 1), verbose)
    , x(x)  // x-coordinate of lower left-hand corner of segment block
    , y(y)  // y-coordinate of lower left-hand corner of segment block
    , z(z)  // horizontal side length of segment block
    , w(w) // vertical side length of block
{}

// Method from SieveBase doesn't give enough primes, so calling the trusty octant sieve.
void SegmentedSieve::setSmallPrimes() {
    if (verbose) {
        cerr << "Calling the OctantSieve to generate smallPrimes..." << endl;
    }
    OctantSieve s((uint64_t)isqrt(maxNorm), false);
    s.run();
    smallPrimes = s.getBigPrimes();
}


void SegmentedSieve::setSieveArray() {
    // sieveArray holds values for Gint's with x <= a < x + z and y <= b < y + w
    if (verbose) {
        cerr << "Building sieve array..." << endl;
    }
    for (uint32_t i = 0; i < z; i++) {
        vector<bool> column(w, true);
        sieveArray.push_back(column);
    }
    if (verbose) {
        printSieveArrayInfo();
    }
}

void SegmentedSieve::crossOffMultiples(gint g) {
    // Let a + bi be the gint and c + di be the co-factor of the multiple we seek.
    // We need the product (a + bi)(c + di) to be contained in the segment block.
    // This gives inequalities x <= ac - bd < x + z and y <= ad + bc < y + w.
    // Solve these for c to get (ax + by) / (a^2 + b^2) <= c <= (a(x+z) + b(y+z)) / (a^2 + b^2).
    // Then d is defined by some max and min conditions.

    // Convert everything to int64_t type
    int64_t a = g.a;
    int64_t b = g.b;
    int64_t N = g.norm();
    int64_t c, cUpper;
    // Trick: to get ceiling of integer division n / m, call (n + m - 1) / m.
    // Trick only works when m, n > 0.
    if (g.b) {
        c = (a * x + b * y + N - 1) / N;  // using trick
        cUpper = (a * (x + z - 1) + b * (y + w - 1)) / N;
    } else {
        c = (x + a - 1) / a;  // using trick
        cUpper = (x + z - 1) / a;
    }
    for (; c <= cUpper; c++) {
        int64_t d, dUpper;  // d isn't necessarily positive
        if (b) {
            // trick doesn't work because of sign; instead using ceil.
            double d1 = double(a * c - x - z + 1) / double(b);
            double d2 = double(y - b * c) / double(a);
            d = int64_t(ceil(max(d1, d2)));
            // Reusing d1 and d2
            d1 = double(a * c - x) / double(b);
            d2 = double(y + w - 1 - b * c) / double(a);
            dUpper = int64_t(floor(min(d1, d2)));
        } else {
            d = (y + a - 1) / a;  // using trick
            dUpper = (y + w - 1) / a;
        }
        int32_t u = a * c - b * d - x;  // u = ac - bd - x
        int32_t v = b * c + a * d - y;  // v = bc + ad - y
        for (; d <= dUpper; d++) {
            sieveArray[u][v] = false;
            u -= b;
            v += a;
        }
    }
}

void SegmentedSieve::setBigPrimes() {
    for (uint32_t a = 0; a <= z; a ++) {
        for (uint32_t b = 0; b <= z; b++) {
            if (sieveArray[a][b]) {
                gint g(a + x, b + y);
                bigPrimes.push_back(g);
            }
        }
    }
}

uint64_t SegmentedSieve::getCountBigPrimes() {
    uint64_t count = 0;
    for (uint32_t a = 0; a <= z; a ++) {
        for (uint32_t b = 0; b <= z; b++) {
            if (sieveArray[a][b]) {
                count++;
            }
        }
    }
    return count;
}