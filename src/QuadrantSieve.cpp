#include <iostream>
#include "../include/QuadrantSieve.hpp"
using namespace std;

void QuadrantSieve::setSmallPrimes() {
    // Trick: Marking all gint's with norm up to sqrt(x) as prime, sorting them, then
    // only using those that are actually prime when we call crossOffMultiples().
    if (verbose) {
        cerr << "Sorting small Gaussian integers by norm..." << endl;
    }
    for (uint32_t a = 1; a <= isqrt(isqrt(maxNorm)); a++) {
        for (uint32_t b = 0; b <= isqrt(isqrt(maxNorm) - a * a); b++) {
            gint g(a, b);
            smallPrimes.push_back(g);
        }
    }
    sort(smallPrimes.begin(), smallPrimes.end());  // sorting by norm (comparison < is built-in to struct gint)
}


void QuadrantSieve::setSieveArray() {
    // Lots of controversy about vector<bool>.
    // Each boolean value is stored as a single bit but pointers are not available.
    // sieveArray holds values for Gint's with a, b >= 0 and a^2 + b^2 <= x.
    if (verbose) {
        cerr << "Building sieve array..." << endl;
    }
    for (uint32_t a = 0; a <= isqrt(x); a++) {
        uint32_t b = isqrt(x - a * a) + 1;
        vector<bool> column(b, true);  // Create a vector of size b with all values true.
        sieveArray.push_back(column);
    }
    sieveArray[0][0] = false;  // 0 is not prime
    sieveArray[1][0] = false;  // 1 is not prime
    sieveArray[0][1] = false;  // i is not prime
    if (verbose) {
        printSieveArrayInfo();
    }
}

void QuadrantSieve::crossOffMultiples(gint g) {
    if (!sieveArray[g.a][g.b]) { return; } // early exit if g isn't actually prime.

    // c + di should range over a full quadrant and satisfy N(c + di)N(a + bi) <= x
    for (uint32_t c = 1; c <= isqrt(x / g.norm()); c++) {  //ignoring c, d = 0, 0
        int32_t u = c * g.a;  // u = ac - bd
        int32_t v = c * g.b;  // v = bc + ad
        for (uint32_t d = 0; d <= isqrt(x / g.norm() - c * c); d++) {
            if (u > 0) {  // u + vi in first quadrant
                sieveArray[u][v] = false;
            } else {  // u + vi in second quadrant; multiplication by -i brings it to first.
                sieveArray[v][-u] = false;
            }
            u -= g.b;
            v += g.a;
        }
    }
    sieveArray[g.a][g.b] = true;  // crossed this off; need to re-mark it as prime
    if (verbose) {
        printProgress(g);
    }
}


void QuadrantSieve::setBigPrimes() {
    if (verbose) {
        cerr << "Gathering primes after sieve..." << endl;
    }
    for (uint32_t a = 1; a <= isqrt(x); a++) {  // Want to stay off imaginary line
        for (uint32_t b = 0; b <= isqrt(x - a * a); b++) {
            if (sieveArray[a][b]) {
                gint g(a, b);
                bigPrimes.push_back(g);
            }
        }
    }
    if (verbose) {
        cerr << "Done with gathering.\n" << endl;
    }
}

uint64_t QuadrantSieve::getCountBigPrimes() {
    if (verbose) {
        cerr << "Counting primes after sieve..." << endl;
    }
    uint64_t count = 0;
    for (uint64_t a = 1; a <= isqrt(x); a++) {  // Want to stay off imaginary line
        for (uint64_t b = 0; b <= isqrt(x - a * a); b++) {
            if (sieveArray[a][b]) {
                count++;
            }
        }
    }
    count *= 4; // four quadrants
    if (verbose) {
        cerr << "Done with count." << endl;
        cerr << "Total number of primes, including associates: " << count << "\n" << endl;
    }
    return count;
}