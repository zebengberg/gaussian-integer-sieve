#include <iostream>
#include "QuadrantSieve.hpp"
using namespace std;

void QuadrantSieve::setSmallPrimes() {
    // Trick: Marking all gint's with norm up to sqrt(x) as prime, sorting them, then
    // only using those that are actually prime when we call crossOffMultiples().
    if (display) {
        cout << "Sorting small Gaussian integers by norm..." << endl;
    }
    for (long a = 1; a <= isqrt(isqrt(maxNorm)); a++) {
        for (long b = 0; b <= isqrt(isqrt(maxNorm) - a * a); b++) {
            gint g(a, b);
            smallPrimes.push_back(g);
        }
    }
    sort(smallPrimes.begin(), smallPrimes.end());  // sorting by norm (this is built-in to struct gint)
}


void QuadrantSieve::setSieveArray() {
    // Lots of controversy about vector<bool>; google it.
    // Each boolean value is stored as a single bit but pointers are not available.
    // sieveArray holds values for Gint's with a, b >= 0 and a^2 + b^2 <= x.
    if (display) {
        cout << "Building sieve array..." << endl;
    }
    for (long a = 0; a <= isqrt(x); a++) {
        long b = isqrt(x - a * a) + 1;
        vector<bool> column((unsigned long)b, true);  // Create a vector of size b with all values true.
        sieveArray.push_back(column);
    }
    sieveArray[0][0] = false;  // 0 is not prime
    sieveArray[1][0] = false;  // 1 is not prime
    sieveArray[0][1] = false;  // i is not prime
    if (display) {
        printSieveArrayInfo();
    }
}

void QuadrantSieve::crossOffMultiples(gint g) {
    if (!sieveArray[g.a][g.b]) { return; } // early exit if g isn't actually prime.
    // c + di should range over a full quadrant and satisfy N(c + di)N(a + bi) <= x
    for (long c = 1; c <= isqrt(x / g.norm()); c++) {  //ignoring c, d = 0, 0
        long u = c * g.a;  // u = ac - bd
        long v = c * g.b;  // v = bc + ad
        for (long d = 0; d <= isqrt(x / g.norm() - c * c); d++) {
            // inexpensive to check this; could avoid this, but would be less understandable
            if (u > 0) {
                sieveArray[u][v] = false;
            } else {  // u + vi in second quadrant; multiplication by -i brings it to first.
                sieveArray[v][-u] = false;
            }
            u -= g.b;
            v += g.a;
        }
    }
    sieveArray[g.a][g.b] = true;  // crossed this off; need to re-mark it as prime
    if (display) {
        printProgress(g);
    }
}

// Old method; keeping for reference.
void QuadrantSieve::crossOffMultiplesAlt(gint g) {
    long p, subgroupSize;
    if (g.b) {  // degree 1 primes
        p = g.norm();
        subgroupSize = p;
    } else {  // degree 2 primes
        p = g.a;
        subgroupSize = 1;  // trivial subgroup
    }
    for (long i = 0; i < subgroupSize; i++) {
        // s + it is an element of the additive subgroup < a + bi > in Z[i] / pZ[i]
        long s = i * g.a % p;
        long t = i * g.b % p;
        for (long u = s; u <= isqrt(x); u += p) {
            for (long v = t; v <= isqrt(x - u * u); v += p) {
                sieveArray[u][v] = false;
            }
        }
    }
    sieveArray[g.a][g.b] = true;  // crossed this off; need to remark it as prime
}

void QuadrantSieve::setBigPrimes() {
    if (display) {
        cout << "Gathering primes after sieve..." << endl;
    }
    for (long a = 1; a <= isqrt(x); a++) {  // Want to stay off imaginary line
        for (long b = 0; b <= isqrt(x - a * a); b++) {
            if (sieveArray[a][b]) {
                gint g(a, b);
                bigPrimes.push_back(g);
            }
        }
    }
}