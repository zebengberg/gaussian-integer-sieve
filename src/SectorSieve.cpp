/* Perform sieving in the sector defined by norm <= x and alpha <= arg < beta.
 * Must have x > 0, alpha >= 0, and beta <= pi/4.
 */



#include <iostream>
#include <cmath>
#include "../include/OctantSieve.hpp"
#include "../include/SectorSieve.hpp"
using namespace std;


// Using an initializer list in the constructor.
SectorSieve::SectorSieve(uint64_t x, double alpha, double beta, bool verbose)
// Calling SieveTemplate constructor to set maxNorm = N(upper right corner)
        : SieveTemplate<bool>(x, verbose)
        , x(x)  // norm bound
        , alpha(alpha)  // beginning angle of sector
        , beta(beta) // final angle of sector
{
    if (beta < alpha) { // the order of alpha and beta doesn't matter; corrected here
        long double temp = alpha;
        this->alpha = beta;  // need this to make reassignment stick outside of this block of code
        this->beta = temp;
        if ((beta > M_PI_2) || (alpha < 0)) {
            cerr << "The interval [alpha, beta) should be a subinterval of [0, pi/2)."
                 << "Choose different values for alpha and beta or use defaults."
                 << endl;
            throw;
        }
    }
}




void SectorSieve::setSmallPrimes() {
    // Setting smallPrimes to the output of QuadrantSieve, which can generate
    // its own small primes.
    if (verbose) {
        cerr << "Calling the OctantSieve to generate smallPrimes..." << endl;
    }
    OctantSieve s((uint64_t)isqrt(maxNorm), false);
    s.run();
    smallPrimes = s.getBigPrimes();
}

// Involves many intersections.
// The annular sector is bounded by the equations
// x^2 + y^2 = maxNorm
// y = tan(alpha) x
// y = tan(beta) x
// Note that the state of the gaussian integer u + iv with in sieveArray is
// given by sieveArray[u][v - floor(tan(alpha) * u)]
void SectorSieve::setSieveArray() {
    // sieveArray holds values for gint's with a, b >= 0
    if (verbose) {
        cerr << "Building sieve array..." << endl;
    }
    // Putting something in at a = 0, and crossing that something off.
    heightShifts.push_back(1);
    sieveArray.emplace_back(1, false);

    for (uint32_t a = 1; a <= isqrt(x / (1 + pow(tan(alpha), 2) )); a++) {
        // Vertical offset from y-axis to sector.
        heightShifts.push_back(int32_t(ceil(tan(alpha) * a - tolerance)));

        // a-value of intersection
        uint32_t intersection = isqrt(x / (1 + pow(tan(beta), 2)));
        int32_t b;  // represents the height of the column at a in sieveArray
        if (a <= intersection) {
            // Do not include any point on the ray theta = beta.
            b = int32_t(floor(tan(beta) * a + tolerance)) - heightShifts[a];
        } else {
            b = isqrt(x - a * a) - int32_t(ceil(tan(alpha) * a - tolerance)) + 1;
        }
        vector<bool> column(b, true);  // Create a vector of size b with all values true.
        sieveArray.push_back(column);
    }
    if (alpha == 0) {
        sieveArray[1][0] = false;  // crossing off 1 (it's not prime)
    }
    if (verbose) {
        printSieveArrayInfo();
    }
}



void SectorSieve::crossOffMultiples(gint g) {
    // First form a mini sector indexed by c and d. Imagining this sector as a
    // slice of pizza, there are three possible configurations for which point
    // on the pizza with tip at origin has largest x-value.
    uint64_t intersectionc1 = isqrt(x / (g.norm() * (1 + pow(tan(beta - g.arg()), 2))));
    uint64_t intersectionc2 = isqrt(x / (g.norm() * (1 + pow(tan(alpha - g.arg()), 2))));
    uint64_t maxIntersectionc = max(intersectionc1, intersectionc2);
    uint64_t cUpper = max(maxIntersectionc, (uint64_t)isqrt(x / g.norm()));
    // Using 64 bits for c since we'll need to square it below.
    for (uint64_t c = 1; c <= cUpper; c++) {  //ignoring c, d = 0, 0
        int32_t d = ceil(max(tan(alpha - g.arg()) * c - tolerance, -(long double)sqrt(x / g.norm() - c * c)));
        // Subtract tolerance since arg(z) strictly less than beta.
        int32_t dUpper = floor(min(tan(beta - g.arg()) * c - tolerance, (long double)sqrt(x / g.norm() - c * c)));
        int32_t u = g.a * c - g.b * d;  // u = ac - bd
        int32_t v = g.b * c + g.a * d;  // v = bc + ad
        for (; d <= dUpper; d++) {
            try {  // round-off error may make this go out of bounds
                sieveArray.at(u).at(uint32_t(v - heightShifts[u])) = false; // get back into sieve array
            } catch (const out_of_range& e) {
                cerr << "Out of range error!" << endl;
                cerr << "a, b = " << g.a << ", " << g.b << endl;
                cerr << "u, v = " << u << ", " << v << endl;
                cerr << "heightShifts[u] = " << heightShifts[u] << endl;
                cerr << "c, cUpper = " << c << ", " << cUpper << endl;
                cerr << "d, dUpper = " << d << ", " << dUpper << endl;
                cerr << "length of column at u: " << sieveArray[u].size() << endl;
                exit(1);
            }
            u -= g.b;
            v += g.a;
        }
    }
    // Checking if gint in sieve array so we can cross it off.
    if ((g.arg() >= alpha) && (g.arg() < beta) && (g.norm() <= x)) {
        sieveArray[g.a][g.b - tan(alpha) * g.a] = true;
    }
    if (verbose) {
        printProgress(g);
    }
}


void SectorSieve::setBigPrimes() {
    if (verbose) {
        cerr << "Gathering primes after sieve..." << endl;
    }
    for (uint32_t a = 0; a <= isqrt(x / (1 + pow(tan(alpha), 2) )); a++) {
        // a-value of intersection
        uint32_t intersection = isqrt(x / (1 + pow(tan(beta), 2)));
        int32_t bUpper;  // represents the height of the column at a in sieveArray
        if (a <= intersection) {
            bUpper = int32_t(tan(beta) * a) - heightShifts[a];
        } else {
            bUpper = isqrt(x - a * a) - heightShifts[a];
        }
        for (int32_t b = 0; b <= bUpper; b++) {
            if (sieveArray[a][b]) {
                gint g(a, uint32_t (b + heightShifts[a]));  // pushing back up into actual sector
                bigPrimes.push_back(g);
            }
        }
    }
    if (verbose) {
        cerr << "Done with gathering.\n" << endl;
    }
}

uint64_t SectorSieve::getCountBigPrimes() {
    if (verbose) {
        cerr << "Counting primes after sieve..." << endl;
    }
    uint64_t count = 0;
    for (uint32_t a = 0; a <= isqrt(x / (1 + pow(tan(alpha), 2) )); a++) {
        // a-value of intersection
        uint32_t intersection = isqrt(x / (1 + pow(tan(beta), 2)));
        int32_t bUpper;  // represents the height of the column at a in sieveArray
        if (a <= intersection) {
            bUpper = int32_t(tan(beta) * a) - int32_t(ceil(tan(alpha) * a));
        } else {
            bUpper = isqrt(x - a * a) - int32_t(ceil(tan(alpha) * a));
        }
        for (int32_t b = 0; b <= bUpper; b++) {
            if (sieveArray[a][b]) {
                count++;
            }
        }
    }
    if (verbose) {
        cerr << "Total number of primes in sector: " << count << "\n" << endl;
    }
    return count;  // four quadrants
}