#include <iostream>
#include <cmath>
#include "../include/OctantSieve.hpp"
#include "../include/QuadrantSieve.hpp"
using namespace std;

void OctantSieve::setSmallPrimes() {
    // Setting smallPrimes to the output of QuadrantSieve, which can generate
    // its own small primes.
    if (verbose) {
        cerr << "Calling the QuadrantSieve to generate smallPrimes..." << endl;
    }
    QuadrantSieve qs((uint64_t)isqrt(maxNorm), false);
    qs.run();
    smallPrimes = qs.getBigPrimes();
}

void OctantSieve::setSieveArray() {
    // sieveArray holds values for gint's with a, b >= 0, a >= b, and a^2 + b^2 <= x.
    if (verbose) {
        cerr << "Building sieve array..." << endl;
    }
    for (uint32_t a = 0; a <= isqrt(x); a++) {
        // Calculating the intersection of circle a^2 + b^2 <= x and the line a = b.
        uint32_t intersection = isqrt(x / 2);
        uint32_t b = a <= intersection ? a + 1 : isqrt(x - a * a) + 1;
        vector<bool> column(b, true);  // Create a vector of size b with all values true.
        sieveArray.push_back(column);
    }
    sieveArray[0][0] = false;  // 0 is not prime
    sieveArray[1][0] = false;  // 1 is not prime
    if (verbose) {
        printSieveArrayInfo();
    }
}

void OctantSieve::crossOffMultiples(gint g) {
    // c + di should range over a full octant and satisfy N(c + di)N(a + bi) <= x
    for (uint32_t c = 1; c <= isqrt(x / g.norm()); c++) {  //ignoring c, d = 0, 0
        int32_t u = c * g.a;  // u = ac - bd
        int32_t v = c * g.b;  // v = bc + ad
        uint32_t intersection = isqrt(x / (2 * g.norm()));
        uint32_t dBound = c <= intersection ? c : isqrt(x / g.norm() - c * c);
        for (long d = 0; d <= dBound; d++) {
            // apply units and conjugate until u + iv is in sieveArray index
            if (u > 0) {
                if (u >= v) {  // u + vi already in first octant
                    sieveArray[u][v] = false;
                } else {  // u + vi in second octant
                    sieveArray[v][u] = false;
                }
            } else {  // u + vi in second quadrant
                if (v >= -u) {  // u + vi in third octant
                    sieveArray[v][-u] = false;
                } else {  // u + vi in fourth octant
                    sieveArray[-u][v] = false;
                }
            }
            u -= g.b;
            v += g.a;
        }
    }
    if (g.a > g.b) {
        sieveArray[g.a][g.b] = true; // crossed this off; need to re-mark it as prime
    } else {
        sieveArray[g.b][g.a] = true; // conjugate crossed this off
    }
    if (verbose) {
        printProgress(g);
    }
}

void OctantSieve::setBigPrimes() {
    if (verbose) {
        cerr << "Gathering primes after sieve..." << endl;
    }
    bigPrimes.emplace_back(1, 1);  // explicitly avoiding ramifying prime 1 + i
    for (uint32_t a = 2; a <= isqrt(x); a++) {
        uint32_t intersection = isqrt(x / 2);
        long bBound = a <= intersection ? a : isqrt(x - a * a);
        for (long b = 0; b <= bBound; b++) {
            if (sieveArray[a][b]) {
                gint g(a, b);
                bigPrimes.push_back(g);
                if (b) {  // prime not on real axis
                    bigPrimes.push_back(g.flip());
                }
            }
        }
    }
    if (verbose) {
        cerr << "Done with gathering.\n" << endl;
    }
}

uint64_t OctantSieve::getCountBigPrimes() {
    if (verbose) {
        cerr << "Counting primes after sieve..." << endl;
    }
    uint64_t count = 1;  // explicitly avoiding ramifying prime 1 + i
    for (uint32_t a = 2; a <= isqrt(x); a++) {
        uint32_t intersection = isqrt(x / 2);
        uint32_t bBound = a <= intersection ? a : isqrt(x - a * a);
        for (long b = 0; b <= bBound; b++) {
            if (sieveArray[a][b]) {
                count++;
                if (b) {  // prime not on real axis
                    count++;
                }
            }
        }
    }
    if (verbose) {
        cerr << "Total number of primes, including associates: " << count << "\n" << endl;
    }
    return 4 * count;  // four quadrants
}