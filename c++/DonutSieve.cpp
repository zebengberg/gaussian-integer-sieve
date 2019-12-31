#include <iostream>
#include <cmath>
#include "DonutSieve.hpp"
#include "SegmentedSieve.hpp"
using namespace std;


void DonutSieve::setMemberVariables() {
    maxNorm = x;
    totalProgress = log(log(maxNorm)) - log(2.0);
}

void DonutSieve::setSmallPrimes() {
    smallPrimes = readPrimesFromFile(isqrt(maxNorm));
    //TODO: remove primes from donut
}

void DonutSieve::setSieveArray() {
    // Every 10 x 10 block of gints should get compressed into a 32-bit structure.
    // We'll use an int for this.
    // Proceed as for octant sieve. The int at (a, b) contains the 10 x 10 block with
    // entries in the grid [10a, 10a + 9] x [10b, 10b + 9].
    cout << "Building sieve array..." << endl;
    for (long a = 0; a <= isqrt(x) / 10; a++) {  // might stick out beyond disk; this is okay
        // Calculating the intersection of circle a^2 + b^2 <= x and the line a = b.
        long intersection = long(sqrt(double(x) / 20.0));
        long b = a <= intersection ? a + 1 : isqrt(x / 100 - a * a) + 1;
        vector<unsigned int> column((unsigned long)b, pow(2, 32) - 1);
        sieveArray.push_back(column);
    }
    // TODO: fix this.
    //sieveArray[0][0] = false;  // 0 is not prime
    //sieveArray[1][0] = false;  // 1 is not prime
}




void DonutSieve::printDonut() {
    SegmentedSieve s(0, 0, 10);
    s.setMemberVariables();
    s.setSieveArray();
    s.crossOffMultiples(gint(1, 1));
    s.crossOffMultiples(gint(2, 1));
    s.crossOffMultiples(gint(1, 2));
    s.printSieveArray();
}