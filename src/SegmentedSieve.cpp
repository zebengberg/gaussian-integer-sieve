#include <iostream>
#include <cmath>
#include "../include/SegmentedSieve.hpp"
#include "../include/OctantSieve.hpp"
using namespace std;

// Using an initializer list
SegmentedSieve::SegmentedSieve(uint32_t x, uint32_t y, uint32_t z, bool verbose)
    : SieveTemplate<bool>((uint64_t)(x + z - 1) * (uint64_t)(y + z - 1), verbose)  // calling SieveTemplate constructor to set maxNorm
    , x(x)  // x-coordinate of lower left-hand corner of segment block
    , y(y)  // y-coordinate of lower left-hand corner of segment block
    , z(z-1)  // side length of segment block; NOTE: originally included x + z, now excluding it.
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
    // sieveArray holds values for Gint's with x <= a <= x + z and y <= b <= y + z
    if (verbose) {
        cerr << "Building sieve array..." << endl;
    }
    for (uint32_t i = 0; i <= z; i++) {
        vector<bool> column(z + 1, true);
        sieveArray.push_back(column);
    }
    if (verbose) {
        printSieveArrayInfo();
    }
}

void SegmentedSieve::crossOffMultiples(gint g) {
    // Let a + bi be the gint and c + di be the co-factor of the multiple we seek.
    // We need the product (a + bi)(c + di) to be contained in the segment block.
    // This gives inequalities x <= ac - bd <= x + z and y <= ad + bc <= y + z.
    // Solve these for c to get (ax + by) / (a^2 + b^2) <= c <= (a(x+z) + b(y+z)) / (a^2 + b^2).
    // Then d is defined by some max and min conditions.

    uint32_t cLower, cUpper;
    // Trick: to get ceiling of integer division n / m, call (n + m - 1) / m.
    // Trick only works when m, n > 0.
    if (g.b) {
        cLower = (uint32_t)(g.a * x + g.b * y + g.norm() - 1) / g.norm();  // using trick
        cUpper = (uint32_t)(g.a * (x + z) + g.b * (y + z)) / g.norm();
    } else {
        cLower = (uint32_t)(x + g.a - 1) / g.a;  // using trick
        cUpper = (uint32_t)(x + z) / g.a;
    }
    for (uint32_t c = cLower; c <= cUpper; c++) {
        int32_t dLower, dUpper;  // d isn't necessarily positive
        if (g.b) {
            // trick doesn't work because of sign; instead using ceil.
            double d1 = double(g.a * c - x - z) / double(g.b);
            double d2 = double(y - g.b * c) / double(g.a);
            dLower = int32_t(ceil(max(d1, d2)));
            d1 = double(g.a * c - x) / double(g.b);
            d2 = double(y + z - g.b * c) / double(g.a);
            dUpper = int32_t(floor(min(d1, d2)));
        } else {
            dLower = (y + g.a - 1) / g.a;  // using trick
            dUpper = (y + z) / g.a;
        }
        int32_t u = g.a * c - g.b * dLower - x;  // u = ac - bd - x
        int32_t v = g.b * c + g.a * dLower - y;  // v = bc + ad - y
        for (int32_t d = dLower; d <= dUpper; d++) {
            sieveArray[u][v] = false;
            u -= g.b;
            v += g.a;
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