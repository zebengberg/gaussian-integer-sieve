#include <iostream>
#include "../include/OctantSieve.hpp"
using namespace std;


// Generate the small primes using the Sieve of Erathosthenes trick in which
// the smallest unmarked integer in the sieve array must be prime. First mark
// all gints in the first octant with norm up to sqrt(x) as prime, sort them,
// then use those that are unmarked calling crossOffMultiples().
void OctantSieve::setSmallPrimes() {
    if (verbose) {
        cerr << "Sorting small Gaussian integers by norm..." << endl;
    }
    // Manually putting in the prime 1 + i.
    smallPrimes.emplace_back(gint(1, 1));
    // Creating a temporary container for small primes.
    vector<gint> tempSmallPrimes;
    // The x-coord of intersection point between the circle
    // a^2 + b^2 <= sqrt(maxNorm) and the line a = b.
    uint32_t intersection = isqrt(isqrt(maxNorm) / 2);
    for (uint32_t a = 2; a <= isqrt(isqrt(maxNorm)); a++) {
        // Not including points on the line a = b.
        uint32_t bUpper = a <= intersection ? a - 1 : isqrt(isqrt(maxNorm) - a * a);
        for (uint32_t b = 0; b <= bUpper; b++) {
            gint g(a, b);
            tempSmallPrimes.push_back(g);
        }
    }
    // Sorting by norm; the comparison < is built-in to struct gint.
    sort(tempSmallPrimes.begin(), tempSmallPrimes.end());
    // Putting the sorted primes and their flipped associates into the real
    // smallPrimes container.
    for (gint g : tempSmallPrimes) {
        smallPrimes.push_back(g);
        if (g.b) {  // Don't want to include imaginary axis.
            smallPrimes.push_back(g.flip());
        }
    }
}


// The sieve array is indexed by gints with b >= 0, a >= b, and
// a^2 + b^2 <= maxNorm.
void OctantSieve::setSieveArray() {
    if (verbose) {
        cerr << "Building sieve array..." << endl;
    }
    // The x-coord of intersection point between circle a^2 + b^2 <= maxNorm
    // and line a = b.
    uint32_t intersection = isqrt(maxNorm / 2);
    for (uint32_t a = 0; a <= isqrt(x); a++) {
        uint32_t height = a <= intersection ? a + 1 : isqrt(maxNorm - a * a) + 1;
        // Create a vector of size height with all values true.
        vector<bool> column(height, true);
        sieveArray.push_back(column);
    }
    sieveArray[0][0] = false;  // 0 is not prime
    sieveArray[1][0] = false;  // 1 is not prime
    if (verbose) {
        printSieveArrayInfo();
    }
}


// If gint g has already been crossed off, then we abort. If not, then g is
// prime, and we cross off its multiples. If g = a + bi, we let c + di be the
// cofactor. So c + di should range over a full octant with norm up to
// maxNorm / g.norm() and satisfy N(c + di) N(a + bi) <= x. Because we are
// working with a full octant, the geometry is easy (as opposed to SectorSieve).
void OctantSieve::crossOffMultiples(gint g) {
    // Early exit if g isn't actually prime.
    if (g.a >= g.b) {
        if (!sieveArray[g.a][g.b]) { return; }
    } else {
        if (!sieveArray[g.b][g.a]) { return; }
    }

    // Ignoring c = 0 and d = 0, which would cross off the gint 0 + 0i.
    for (uint32_t c = 1; c <= isqrt(maxNorm / g.norm()); c++) {
        int32_t u = c * g.a;  // u = ac - bd
        int32_t v = c * g.b;  // v = bc + ad
        uint32_t intersection = isqrt(maxNorm / (2 * g.norm()));
        uint32_t dUpper = c <= intersection ? c : isqrt(maxNorm / g.norm() - c * c);
        for (uint32_t d = 0; d <= dUpper; d++) {
            // Apply units and conjugate until u + vi is in index of sieveArray.
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
    // Uncrossing the gint itself and its conjugate which we crossed off when
    // c = 1 and d = 0.
    if (g.a > g.b) {
        sieveArray[g.a][g.b] = true;
    } else {
        sieveArray[g.b][g.a] = true;
    }
    if (verbose) {
        printProgress(g);
    }
}


// Stepping through sieveArray and storing unmarked entries as gints.
void OctantSieve::setBigPrimes() {
    if (verbose) {
        cerr << "Gathering primes after sieve..." << endl;
    }
    // Explicitly avoiding ramifying prime 1 + i
    if (maxNorm >= 2) { bigPrimes.emplace_back(1, 1); }

    uint32_t intersection = isqrt(maxNorm / 2);
    for (uint32_t a = 2; a <= isqrt(maxNorm); a++) {
        // Can avoid line a = b.
        uint32_t bUpper = a <= intersection ? a - 1 : isqrt(maxNorm - a * a);
        for (uint32_t b = 0; b <= bUpper; b++) {
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


// Counting primes after sieve and returning the count.
uint64_t OctantSieve::getCountBigPrimes() {
    if (verbose) {
        cerr << "Counting primes after sieve..." << endl;
    }
    if (maxNorm < 2) { return 0; }
    uint64_t count = 1;  // manually count 1 + i

    uint32_t intersection = isqrt(maxNorm / 2);
    for (uint32_t a = 2; a <= isqrt(maxNorm); a++) {
        // Can avoid line a = b.
        uint32_t bUpper = a <= intersection ? a - 1 : isqrt(maxNorm - a * a);
        for (uint32_t b = 0; b <= bUpper; b++) {
            if (sieveArray[a][b]) {
                count++;
                if (b) {  // prime not on real axis; count the flipped version
                    count++;
                }
            }
        }
    }
    count *= 4;  // four quadrants
    if (verbose) {
        cerr << "Total number of primes, including associates: " << count << "\n" << endl;
    }
    return count;
}