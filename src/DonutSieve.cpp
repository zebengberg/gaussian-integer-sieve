#include <iostream>
#include <cmath>
#include "../include/DonutSieve.hpp"
#include "../include/QuadrantSieve.hpp"
using namespace std;


DonutSieve::DonutSieve(long x, bool display)
    : SieveTemplate<unsigned int>(x, display)  // calling SieveTemplate constructor to set maxNorm
    , x(x)

    // For a fixed c, d will iterate by d += 2 or d += 4, similar to the classic
    // wheel sieve. Got these iteration patterns in printDonut(), then pasted them
    // here so they can be hard coded in constructor.

    // dStart array determines the starting value of d for a given c during crossOffMultiples()
    , dStart{1, 0, 3, 0, 1, 2, 1, 0, 3, 0}

    // donut array tracks how to move from one d to another during crossOffMultiples()
    // note: rows and columns reversed from usual complex plane orientation -- don't get confused.
    , gapDonut{{0, 2, 0, 4, 0, 0, 0, 2, 0, 2},
               {4, 0, 0, 0, 2, 0, 4, 0, 0, 0},
               {0, 0, 0, 2, 0, 2, 0, 6, 0, 0},
               {2, 0, 6, 0, 0, 0, 0, 0, 2, 0},
               {0, 4, 0, 0, 0, 4, 0, 0, 0, 2},
               {0, 0, 2, 0, 2, 0, 2, 0, 4, 0},
               {0, 4, 0, 0, 0, 4, 0, 0, 0, 2},
               {2, 0, 6, 0, 0, 0, 0, 0, 2, 0},
               {0, 0, 0, 2, 0, 2, 0, 6, 0, 0},
               {4, 0, 0, 0, 2, 0, 4, 0, 0, 0}}

    , bitDonut{{99, 0, 99, 1, 99, 99, 99, 2, 99, 3},
               {4, 99, 99, 99, 5, 99, 6, 99, 99, 99},
               {99, 99, 99, 7, 99, 8, 99, 9, 99, 99},
               {10, 99, 11, 99, 99, 99, 99, 99, 12, 99},
               {99, 13, 99, 99, 99, 14, 99, 99, 99, 15},
               {99, 99, 16, 99, 17, 99, 18, 99, 19, 99},
               {99, 20, 99, 99, 99, 21, 99, 99, 99, 22},
               {23, 99, 24, 99, 99, 99, 99, 99, 25, 99},
               {99, 99, 99, 26, 99, 27, 99, 28, 99, 99},
               {29, 99, 99, 99, 30, 99, 31, 99, 99, 99}}

    , realPartDecompress{0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9}
    , imagPartDecompress{1, 3, 7, 9, 0, 4, 6, 3, 5, 7, 0, 2, 8, 1, 5, 9, 2, 4, 6, 8, 1, 5, 9, 0, 2, 8, 3, 5, 7, 0, 4, 6}
{}

void DonutSieve::setSmallPrimes() {
    // Setting smallPrimes to the output of QuadrantSieve, which can generate
    // its own small primes.
    if (verbose) {
        cerr << "Calling the QuadrantSieve to generate smallPrimes..." << endl;
    }
    QuadrantSieve qs(isqrt(maxNorm), false);
    qs.run();
    smallPrimes = qs.getBigPrimes();
}

void DonutSieve::setSieveArray() {
    // Every 10 x 10 block of gints should get compressed into a 32-bit structure.
    // We'll use an int for this.
    // Proceed as for octant sieve. The int at (a, b) contains the 10 x 10 block with
    // entries in the grid [10a, 10a + 9] x [10b, 10b + 9].

    // Sieve array might stick out beyond boundary of disk -- we'll fix this in getBigPrimes.
    if (verbose) {
        cerr << "Building sieve array..." << endl;
    }
    for (long a = 0; a <= isqrt(x) / 10; a++) {
        // Calculating the intersection of circle a^2 + b^2 <= x and the line a = b.
        long intersection = long(sqrt(double(x) / 200.0));
        long b = a <= intersection ? a + 1: isqrt(x / 100 - a * a) + 1;
        vector<unsigned int> column((unsigned long)b, (unsigned int)pow(2, 32) - 1);
        sieveArray.push_back(column);
    }
    setFalse(1, 0);  // 1 is not prime
    setFalse(0, 1);  // i is not prime
    if (verbose) {
        printSieveArrayInfo();
    }
}

void DonutSieve::crossOffMultiples(gint g) {
    if (g.norm() <= 5) { return; }  // exit early if g is above 2 or 5.
    // Let a + bi be the gint and c + di be the co-factor of the multiple we seek.
    // Because the product (a + bi)(c + di) should be coprime to 10, we need that
    // c + di is also coprime to 10. This gives conditions on c and d mod 10.

    // As with all of these crossOffMultiples() methods, we use a double for-loop
    // of the form for (iterate over c's) { for (iterate over d's) }.
    for (long c = 0; c <= isqrt(x / g.norm()); c++) {
        long d = dStart[c % 10];  // starting value of d for while loop
        long u = c * g.a - d * g.b;  // u = ac - bd
        long v = c * g.b + d * g.a;  // v = bc + ad

        long intersection = long(sqrt(double(x) / double(2 * g.norm())));
        long dBound = c <= intersection ? c : isqrt(x / g.norm() - c * c);
        long jump;
        while (d <= dBound) {  // replacing inner for loop with while loop
            // apply units and conjugate until u + iv is in sieveArray index
            if (u > 0) {
                if (u >= v) {  // u + vi already in first octant
                    setFalse(u, v);
                } else {  // u + vi in second octant
                    setFalse(v, u);
                }
            } else {  // u + vi in second quadrant
                if (v >= -u) {  // u + vi in third octant
                    setFalse(v, -u);
                } else {  // u + vi in fourth octant
                    setFalse(-u, v);
                }
            }
            jump = gapDonut[c % 10][d % 10];
            u -= jump * g.b;
            v += jump * g.a;
            d += jump;
        }
    }
    if (g.a > g.b) {
        setTrue(g.a, g.b); // crossed this off; need to re-mark it as prime
    } else {
        setTrue(g.b, g.a);
    }
}

void DonutSieve::setFalse(long u, long v) {
    // Set the correct bit in the sieveArray to false corresponding to the gint u + vi
    unsigned int bit = bitDonut[u % 10][v % 10];
    sieveArray[u / 10][v / 10] &= ~(1u << bit);  // clearing the bit; 1u is unsigned int
}

void DonutSieve::setTrue(long u, long v) {
    // Set the correct bit in the sieveArray to false corresponding to the gint u + vi
    unsigned int bit = bitDonut[u % 10][v % 10];
    sieveArray[u / 10][v / 10] |= (1u << bit);  // setting the bit to 0; 1u is unsigned int
}

void DonutSieve::setBigPrimes() {
    if (verbose) {
        cerr << "Gathering primes after sieve..." << endl;
    }
    // Putting in primes dividing 10.
    bigPrimes.emplace_back(1, 1);
    bigPrimes.emplace_back(2, 1);
    bigPrimes.emplace_back(1, 2);
    for (long a = 0; a <= isqrt(x) / 10; a++) {
        long intersection = long(sqrt(double(x) / 20.0));
        long bBound = a <= intersection ? a : isqrt(x / 100 - a * a);
        for (long b = 0; b <= bBound; b++) {
            for (unsigned int bit = 0; bit < 32; bit++) {
                if ((sieveArray[a][b] >> bit) & 1u) {
                    gint g(10 * a + realPartDecompress[bit], 10 * b + imagPartDecompress[bit]);
                    // check for boundary blocks and to avoid imag multiple of degree 2
                    if ((g.norm() <= x) && (g.a) && (g.a > g.b)) {
                        bigPrimes.push_back(g);
                        if (g.b) {  // prime not on real axis
                            bigPrimes.push_back(g.flip());
                        }
                    }
                }
            }
        }
    }
    if (verbose) {
        cerr << "Done with gathering.\n" << endl;
    }
}

unsigned long DonutSieve::getCountBigPrimes() {
    if (verbose) {
        cerr << "Counting primes after sieve..." << endl;
    }
    unsigned long count = 3;  // 3 primes dividing 10
    for (long a = 0; a <= isqrt(x) / 10; a++) {
        long intersection = long(sqrt(double(x) / 20.0));
        long bBound = a <= intersection ? a : isqrt(x / 100 - a * a);
        for (long b = 0; b <= bBound; b++) {
            for (unsigned int bit = 0; bit < 32; bit++) {
                if ((sieveArray[a][b] >> bit) & 1u) {
                    // Coordinates of actual gint.
                    long aa = 10 * a + realPartDecompress[bit];
                    long bb = 10 * b + imagPartDecompress[bit];
                    // check for boundary blocks and to avoid imag multiple of degree 2
                    if ((aa * aa + bb * bb <= x) && aa && (aa > bb)) {
                        count++;
                        if (bb) {  // prime not on real axis
                            count++;
                        }
                    }
                }
            }
        }
    }
    if (verbose) {
        cerr << "Done with count." << endl;
        cerr << "Total number of primes, including associates: " << count << "\n" << endl;
    }
    return 4 * count;  // four quadrants
}



// Using this method to print out some of the member variables used in constructor.
// This prints out necessary data for the donut.
void DonutSieve::printDonutArrays() {
    SegmentedSieve s(0, 0, 15, false);  // going a little bit beyond 9 so we can get gaps
    s.setSieveArray();
    s.crossOffMultiples(gint(1, 1));
    s.crossOffMultiples(gint(2, 1));
    s.crossOffMultiples(gint(1, 2));

    // Printing the gapDonut. This array contains information about the horizontal space of the donut
    // and is used to iterate d when calling crossOffMultiples().
    cout << "gapDonut:" << endl;
    for (unsigned int c = 0; c < 10; c++) {
        string line("{");
        for (unsigned int d = 0; d < 10; d++) {
            if (s.getSieveArrayValue(c, d)) {
                unsigned int e = 0;
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

    // Printing the dStart array. For a given c, this array gives the starting value of d.
    cout << "\n\ndStart:" << endl;
    string arr = "{";
    for (unsigned int c = 0; c < 10; c++) {
        unsigned int d = 0;
        while (!s.getSieveArrayValue(c, d)) { d++; }
        arr += to_string(d);
        arr += ", ";
    }
    arr.pop_back();
    arr.pop_back();
    arr += "}";
    cout << arr << endl;

    // Printing the bitDonut array. Used to compress a residue mod 10Z[i] into a bit position.
    cout << "\n\nbitDonut:" << endl;
    int bit = 0;
    for (unsigned int c = 0; c < 10; c++) {
        string line = "{";
        for (unsigned int d = 0; d < 10; d++) {
            if (s.getSieveArrayValue(c, d)) {
                line += to_string(bit);
                line += ", ";
                bit++;
            } else {
                line += "99, ";  // put in some crazy value that'll hopefully throw an error if accessed
            }
        }
        line.pop_back();
        line.pop_back();
        line += "},";
        cout << line << endl;
    }

    // Printing the decompression array, which is used in setBigPrimes().
    unsigned int c = 0;
    unsigned int d = 0;
    string real = "{";
    string imag = "{";
    for (bit = 0; bit < 32; bit++) {
        do {
            if (d < 9) {
                d++;
            } else {
                c++;
                d = 0;
            }
        } while (!s.getSieveArrayValue(c, d));
        real += to_string(c);
        real += ", ";
        imag += to_string(d);
        imag += ", ";
    }
    cout << "\n\nrealPartDecompress:" << endl;
    real.pop_back();
    real.pop_back();
    real += "}";
    cout << real << endl;
    cout << "imagPartDecompress:" << endl;
    imag.pop_back();
    imag.pop_back();
    imag += "}";
    cout << imag << endl;
}