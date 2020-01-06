#include <iostream>
#include <cmath>
#include "SegmentedSieve.hpp"
#include "OctantSieve.hpp"
using namespace std;

// Using an initializer list
SegmentedSieve::SegmentedSieve(long x, long y, long z, bool display)
    : SieveTemplate<bool>((x + z) * (y + z), display)  // calling SieveTemplate constructor to set maxNorm
    , x(x)  // x-coordinate of lower left-hand corner of segment block
    , y(y)  // y-coordinate of lower left-hand corner of segment block
    , z(z)  // side length of segment block
{}

// Method from SieveBase doesn't give enough primes, so calling the trusty octant sieve.
void SegmentedSieve::setSmallPrimes() {
    if (display) {
        cout << "Calling the OctantSieve to generate smallPrimes..." << endl;
    }
    OctantSieve s(isqrt(maxNorm), false);
    smallPrimes = s.run();
}


void SegmentedSieve::setSieveArray() {
    // sieveArray holds values for Gint's with x <= a <= x + z and y <= b <= y + z
    if (display) {
        cout << "Building sieve array..." << endl;
    }
    for (long i = 0; i <= z; i ++) {
        vector<bool> column((unsigned long)(z + 1), true);
        sieveArray.push_back(column);
    }
    if (display) {
        printSieveArrayInfo();
    }
}

void SegmentedSieve::crossOffMultiples(gint g) {
    // Let a + bi be the gint and c + di be the co-factor of the multiple we seek.
    // We need the product (a + bi)(c + di) to be contained in the segment block.
    // This gives inequalities x <= ac - bd <= x + z and y <= ad + bc <= y + z.
    // Solve these for c to get (ax + by) / (a^2 + b^2) <= c <= (a(x+z) + b(y+z)) / (a^2 + b^2).
    // Then d is defined by some max and min conditions.

    long cLower, cUpper;
    // Trick: to get ceiling of integer division n / m, call (n + m - 1) / m.
    // Trick only works when m, n > 0.
    if (g.b) {
        cLower = (g.a * x + g.b * y + g.norm() - 1) / g.norm();  // using trick
        cUpper = (g.a * (x + z) + g.b * (y + z)) / g.norm();
    } else {
        cLower = (x + g.a - 1) / g.a;  // using trick
        cUpper = (x + z) / g.a;
    }
    for (long c = cLower; c <= cUpper; c++) {
        long dLower, dUpper;
        if (g.b) {
            // trick doesn't work because of sign; instead using ceil.
            double d1 = double(g.a * c - x - z) / double(g.b);
            double d2 = double(y - g.b * c) / double(g.a);
            dLower = long(ceil(max(d1, d2)));
            d1 = double(g.a * c - x) / double(g.b);
            d2 = double(y + z - g.b * c) / double(g.a);
            dUpper = long(floor(min(d1, d2)));
        } else {
            dLower = (y + g.a - 1) / g.a;  // using trick
            dUpper = (y + z) / g.a;
        }
        long u = g.a * c - g.b * dLower - x;  // u = ac - bd - x
        long v = g.b * c + g.a * dLower - y;  // v = bc + ad - y
        for (long d = dLower; d <= dUpper; d++) {
            sieveArray[u][v] = false;
            u -= g.b;
            v += g.a;
        }
    }
}

void SegmentedSieve::setBigPrimes() {
    for (long a = 0; a <= z; a ++) {
        for (long b = 0; b <= z; b++) {
            if (sieveArray[a][b]) {
                gint g(a + x, b + y);
                bigPrimes.push_back(g);
            }
        }
    }
}