#include <iostream>
#include <cmath>
#include "BlockSieve.hpp"
#include "../include/OctantSieve.hpp"
using namespace std;

// Using an initializer list
BlockSieve::BlockSieve(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy, bool verbose)
    // Calling SieveTemplate constructor to set maxNorm
    : SieveTemplate<bool>((uint64_t)(x + dx - 1) * (uint64_t)(y + dy - 1), verbose)
    , x(x)  // x-coordinate of lower left-hand corner of segment block
    , y(y)  // y-coordinate of lower left-hand corner of segment block
    , dx(dx)  // horizontal side length of segment block
    , dy(dy) // vertical side length of block
{}

// Method from SieveBase doesn't give enough primes, so calling the trusty octant sieve.
void BlockSieve::setSmallPrimes() {
    if (verbose) {
        cerr << "Calling the OctantSieve to generate smallPrimes..." << endl;
    }
    OctantSieve s((uint64_t)isqrt(maxNorm), false);
    s.run();
    smallPrimes = s.getBigPrimes();
}

// Sieve array holds indices corresponding to gints with x <= a < x + dx and
// y <= b < y + dy. The size of the sieve array is dx by dy.
void BlockSieve::setSieveArray() {
    if (verbose) {
        cerr << "Building sieve array..." << endl;
    }
    for (uint32_t i = 0; i < dx; i++) {
        vector<bool> column(dy, true);
        sieveArray.push_back(column);
    }
    if ((x <= 1) && (y == 0)) {
        sieveArray[1 - x][0] = false;  // Crossing off 1
    }
    if ((y <= 1) && (x == 0)) {
        sieveArray[0][1 - y] = false;  // Crossing off i
    }
    if (verbose) {
        printSieveArrayInfo();
    }
}

// Cross off multiples of the gint g = a + bi within the sieveArray.
// Let a + bi be the gint and c + di be the co-factor of the multiple we seek.
// We must over integers c and d such that the product (a + bi)(c + di) is
// contained in the segment block. This gives inequalities x <= ac - bd < x + dx
// and y <= ad + bc < y + dy. Solve these for c to get
// (ax + by) / (a^2 + b^2) <= c <= (a(x + dx - 1) + b(y + dy - 1)) / (a^2 + b^2).
// Then d is defined by some max and min conditions.
void BlockSieve::crossOffMultiples(gint g) {

    // Convert everything to int64_t type
    int64_t a = g.a;
    int64_t b = g.b;
    int64_t N = g.norm();
    int64_t c, cUpper;
    // Trick: to get ceiling of integer division n / m, call (n + m - 1) / m.
    // This only works when m, n > 0.
    if (g.b) {
        c = (a * x + b * y + N - 1) / N;  // using trick
        cUpper = (a * (x + dx - 1) + b * (y + dy - 1)) / N;
    } else {
        c = (x + a - 1) / a;  // using trick
        cUpper = (x + dx - 1) / a;
    }
    for (; c <= cUpper; c++) {
        int64_t d, dUpper;  // d isn't necessarily positive
        if (b) {
            // trick doesn't work because of sign; instead using ceil.
            double d1 = double(a * c - x - dx + 1) / double(b);
            double d2 = double(y - b * c) / double(a);
            d = int64_t(ceil(max(d1, d2)));
            // Reusing d1 and d2
            d1 = double(a * c - x) / double(b);
            d2 = double(y + dy - 1 - b * c) / double(a);
            dUpper = int64_t(floor(min(d1, d2)));
        } else {
            d = (y + a - 1) / a;  // using trick
            dUpper = (y + dy - 1) / a;
        }
        int32_t u = a * c - b * d - x;  // u = ac - bd - x
        int32_t v = b * c + a * d - y;  // v = bc + ad - y
        for (; d <= dUpper; d++) {
            sieveArray[u][v] = false;
            u -= b;
            v += a;
        }
    }
    if ((x <= a) && (a < x + dx) && (y <= b) && (b < y + dy)){
        // crossed this off; need to re-mark it as prime
        sieveArray[a - x][b - y] = true;
    }
    if ((x <= b) && (b < x + dx) && (y <= a) && (a < y + dy)){
        // crossed this off; need to re-mark it as prime
        sieveArray[b - x][a - y] = true;
    }
}

// Combing through the sieve array to find primes after sieving.
void BlockSieve::setBigPrimes() {
    if (verbose) {
        cerr << "Gathering primes after sieve..." << endl;
    }
    for (uint32_t a = 0; a < dx; a ++) {
        for (uint32_t b = 0; b < dy; b++) {
            if (sieveArray[a][b]) {
                gint g(a + x, b + y);
                bigPrimes.push_back(g);
            }
        }
    }
    if (verbose) {
        cerr << "Done gathering." << endl;
    }
}

uint64_t BlockSieve::getCountBigPrimes() {
    if (verbose) {
        cerr << "Counting primes after sieve..." << endl;
    }
    uint64_t count = 0;
    for (uint32_t a = 0; a < dx; a ++) {
        for (uint32_t b = 0; b < dy; b++) {
            if (sieveArray[a][b]) {
                count++;
            }
        }
    }
    if (verbose) {
        cerr << "Total number of primes: " << count << "\n" << endl;
    }
    return count;
}