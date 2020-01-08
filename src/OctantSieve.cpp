#include <iostream>
#include <cmath>
#include "OctantSieve.hpp"
#include "QuadrantSieve.hpp"
using namespace std;

void OctantSieve::setSmallPrimes() {
    // Setting smallPrimes to the output of QuadrantSieve, which can generate
    // its own small primes.
    if (display) {
        cout << "Calling the QuadrantSieve to generate smallPrimes..." << endl;
    }
    QuadrantSieve qs(isqrt(maxNorm), false);
    qs.run();
    smallPrimes = qs.getBigPrimes();
}

void OctantSieve::setSieveArray() {
    // sieveArray holds values for Gint's with a, b >= 0, a >= b, and a^2 + b^2 <= x.
    if (display) {
        cout << "Building sieve array..." << endl;
    }
    for (long a = 0; a <= isqrt(x); a++) {
        // Calculating the intersection of circle a^2 + b^2 <= x and the line a = b.
        long intersection = long(sqrt(double(x) / 2.0));
        long b = a <= intersection ? a + 1 : isqrt(x - a * a) + 1;
        vector<bool> column((unsigned long)b, true);  // Create a vector of size b with all values true.
        sieveArray.push_back(column);
    }
    sieveArray[0][0] = false;  // 0 is not prime
    sieveArray[1][0] = false;  // 1 is not prime
    if (display) {
        printSieveArrayInfo();
    }
}

void OctantSieve::crossOffMultiples(gint g) {
    // c + di should range over a full octant and satisfy N(c + di)N(a + bi) <= x
    for (long c = 1; c <= isqrt(x / g.norm()); c++) {  //ignoring c, d = 0, 0
        long u = c * g.a;  // u = ac - bd
        long v = c * g.b;  // v = bc + ad
        long intersection = long(sqrt(double(x) / double(2 * g.norm())));
        long dBound = c <= intersection ? c : isqrt(x / g.norm() - c * c);
        for (long d = 0; d <= dBound; d++) {
            // apply units and conjugate until u + iv is in sieveArray index
            if (u > 0) {
                if (u >= v) {
                    sieveArray[u][v] = false;  // u + vi already in first octant
                } else {
                    sieveArray[v][u] = false;  // u + vi in second octant
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
    if (display) {
        printProgress(g);
    }
}

void OctantSieve::setBigPrimes() {
    if (display) {
        cout << "Gathering primes after sieve..." << endl;
    }
    bigPrimes.emplace_back(1, 1);  // explicitly avoiding ramifying prime 1 + i
    for (long a = 2; a <= isqrt(x); a++) {
        long intersection = long(sqrt(double(x) / 2.0));
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
    if (display) {
        cout << "Done with gathering.\n" << endl;
    }
}

unsigned long OctantSieve::getCountBigPrimes() {
    if (display) {
        cout << "Counting primes after sieve..." << endl;
    }
    unsigned long count = 1;  // explicitly avoiding ramifying prime 1 + i
    for (long a = 2; a <= isqrt(x); a++) {
        long intersection = long(sqrt(double(x) / 2.0));
        long bBound = a <= intersection ? a : isqrt(x - a * a);
        for (long b = 0; b <= bBound; b++) {
            if (sieveArray[a][b]) {
                count++;
                if (b) {  // prime not on real axis
                    count++;
                }
            }
        }
    }
    if (display) {
        cout << "Done with count." << endl;
        cout << "Total number of primes, including associates: " << count << "\n" << endl;
    }
    return 4 * count;  // four quadrants
}