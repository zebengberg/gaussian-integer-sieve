#include <iostream>
#include <cmath>
#include "DonutSieve.hpp"
#include "SegmentedSieve.hpp"
using namespace std;


DonutSieve::DonutSieve(long x)
    : x(x)
    , SieveTemplate<unsigned int>(x)  // calling SieveTemplate constructor to set maxNorm

    // For a fixed c, d will iterate by d += 2 or d += 4, similar to the classic
    // wheel sieve. Got these iteration patterns in printDonut(), then pasted them
    // here so they can be hard coded in constructor.

    // dStart array determines the starting value of d for a given c during crossOffMultiples()
    , dStart{1, 0, 3, 0, 1, 2, 1, 0, 3, 0}

    // donut array tracks how to move from one d to another during crossOffMultiples()
    // note: rows and columns reversed from usual complex plane orientation -- don't get confused.
    , donut{{0, 2, 0, 4, 0, 0, 0, 2, 0, 2},
            {4, 0, 0, 0, 2, 0, 4, 0, 0, 0},
            {0, 0, 0, 2, 0, 2, 0, 6, 0, 0},
            {2, 0, 6, 0, 0, 0, 0, 0, 2, 0},
            {0, 4, 0, 0, 0, 4, 0, 0, 0, 2},
            {0, 0, 2, 0, 2, 0, 2, 0, 4, 0},
            {0, 4, 0, 0, 0, 4, 0, 0, 0, 2},
            {2, 0, 6, 0, 0, 0, 0, 0, 2, 0},
            {0, 0, 0, 2, 0, 2, 0, 6, 0, 0},
            {4, 0, 0, 0, 2, 0, 4, 0, 0, 0}}
{}

void DonutSieve::setSmallPrimes() {
    smallPrimes = readPrimesFromFile(isqrt(maxNorm));
    //TODO: remove primes above 2 and 5
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
    SegmentedSieve s(0, 0, 15);  // going a little bit beyond 9 so we can get gaps
    s.setSieveArray();
    s.crossOffMultiples(gint(1, 1));
    s.crossOffMultiples(gint(2, 1));
    s.crossOffMultiples(gint(1, 2));

    // Printing the donut and the horizontal gap sizes inside the donut
    for (int c = 0; c < 10; c++) {
        string line("{");
        for (int d = 0; d < 10; d++) {
            if (s.getSieveArrayValue(c, d)) {
                int e = 0;
                do { e++; } while(!s.getSieveArrayValue(c, d + e));
                line += to_string(e);
                line += ", ";
            } else {
                line += "0, ";
            }
        }
        line.pop_back();
        line.pop_back();
        line += "},";
        cout << line << endl;
    }

    // For each c, getting the starting value of d.
    string arr("{");
    for (int c = 0; c < 10; c++) {
        int d = 0;
        while (!s.getSieveArrayValue(c, d)) { d++; }
        arr += to_string(d);
        arr += ", ";
    }
    arr.pop_back();
    arr.pop_back();
    arr += "}";
    cout << "\n\n" << arr << endl;
}

void DonutSieve::crossOffMultiples(gint g) {
    // Let a + bi be the gint and c + di be the co-factor of the multiple we seek.
    // Because the product (a + bi)(c + di) should be coprime to 10, we need that
    // c + di is also coprime to 10. This gives conditions on c and d mod 10.

    // As with all of these crossOffMultiples() methods, we use a double for-loop
    // of the form for (iterate over c's) { for (iterate over d's) }.


    for (long c = 1; c <= isqrt(x / g.norm()); c++) {  // ignoring c, d = 0, 0

        long d = dStart[c % 10];  // starting value of d for while loop

        long u = c * g.a - d * g.b;  // u = ac - bd
        long v = c * g.b + d * g.a;  // v = bc + ad
        long intersection = long(sqrt(double(x) / double(2 * g.norm())));
        long dBound = c <= intersection ? c : isqrt(x / g.norm() - c * c);
        while (d <= dBound) {  // replacing inner for loop with while loop
            d += donut[c % 10][d % 10];
        }
    }


}

void DonutSieve::setBigPrimes() {}