#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "../include/BaseSieve.hpp"

// Will call this constructor from derived classes.
SieveBase::SieveBase(uint64_t maxNorm, bool verbose)
    // initializer list
    : maxNorm(maxNorm)
    , verbose(verbose)
    {
    progress = 0.0;
    // Using PNT to make approximate progress bar.
    totalProgress = log(log(maxNorm)) - log(2.0);
}

void SieveBase::sieve() {
    if (verbose) {
        cerr << "Starting to sieve..." << endl;
    }
    auto startTime = chrono::high_resolution_clock::now();
    for (gint g : smallPrimes) {
        crossOffMultiples(g);
    }
    auto endTime = chrono::high_resolution_clock::now();
    auto totalTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);
    double printTime = double(totalTime.count()) / 1000.0;
    if (verbose) {
        cerr << "Done sieving. Total time for sieving: " << printTime << " seconds." << endl;
    }
}

void SieveBase::printProgress(gint g) {
    int barSize = 80;
    progress += 1.0 / double(g.norm());
    int barPos = int(double(barSize) * progress / totalProgress);
    cerr << "[";
    for (int i = 0; i < barPos; i++) { cerr << "|"; }
    for (int i = barPos; i <= barSize; i++) { cerr << " "; }
    cerr << "]\r";
    cerr.flush();
}

void SieveBase::run() {
    setSmallPrimes();
    setSieveArray();
    sieve();
}

vector<gint> SieveBase::getBigPrimes() {
    setBigPrimes();
    sortBigPrimes();
    return bigPrimes;
}

void SieveBase::sortBigPrimes() {
    if (verbose) {
        cerr << "Sorting primes by norm..." << endl;
    }
    sort(bigPrimes.begin(), bigPrimes.end()); }

void SieveBase::printBigPrimes() {
    for (gint g : bigPrimes) {
        cout << g.a << " " << g.b << endl;
    }
    cerr << "Total number of primes, including associates: " << 4 * bigPrimes.size() << endl;
}

void SieveBase::writeBigPrimesToFile() {
    ofstream f;
    f.open("cpp_primes.csv");
    for (gint g : bigPrimes) {
        f << g.a << ' ' << g.b << endl;
    }
    f.close();
}

// Getting all primes from file small_primes.txt with norm up to maxNorm variable.
// Old method; keeping as a reference.
void SieveBase::setSmallPrimesFromFile() {
    // If maxNorm is too small, then we are not actually taking any primes.
    if (maxNorm < 2) {
        cerr << "No primes read and no sieving needed.";
        exit(1);
    }
    // small_primes.txt contains a list of Gaussian primes sorted by norm.
    ifstream f;
    f.open("PUT IN PATH TO FILE");
    if (!f) {
        cerr << "Unable to open file.";
        exit(1);
    } else {
        cerr << "Reading small primes from file." << endl;
    }

    int64_t a, b;
    f >> a >> b;  // token-based parsing
    gint g(a, b);
    // Need primes up to square root of the maxNorm
    while (g.norm() <= isqrt(maxNorm)) {
        smallPrimes.push_back(g);
        f >> a >> b;  // reading next pair a, b from file f
        g = gint(a, b);
        // If we get to the end of the file, then we don't have enough precomputed primes.
        if (f.eof()) {
            cerr << "Not enough primes in PATH TO FILE" << endl;
            cerr << "Need to override this method with a sieve instance." << endl;
            f.close();
            exit(1);
        }
    }
    f.close();
}


// Specializations of SieveTemplate methods

template <>
void SieveTemplate<bool>::printSieveArrayInfo() {
    uint64_t totalSize = sizeof(sieveArray);
    for (const auto& column : sieveArray) {
        totalSize += sizeof(column) + column.capacity() / 8;  // each bool stored as a bit
    }
    totalSize /= pow(10, 6);  // convert to MB
    cerr << "Sieve array approximate memory use: " << totalSize  << "MB." << endl;
}

template <>
void SieveTemplate<uint32_t>::printSieveArrayInfo() {
    uint64_t totalSize = sizeof(sieveArray);
    for (const auto& column : sieveArray) {
        totalSize += sizeof(column) + column.capacity() * sizeof(unsigned int);
    }
    totalSize /= pow(10, 6);  // convert to MB
    cerr << "Sieve array approximate memory use: " << totalSize  << "MB." << endl;
}

template <>
void SieveTemplate<bool>::printSieveArray() {
    // Print sieve array with same orientation as that in the complex plane.
    uint32_t columnMaxSize = 0;
    for (auto &column : sieveArray) {
        if (column.size() > columnMaxSize) {
            columnMaxSize = column.size();
        }
    }
    // Some type casting because subtraction confuses it.
    for (auto v = (int32_t)columnMaxSize - 1; v >= 0; v--) {
        string row;
        for (auto &column : sieveArray) {
            if (column.size() > v) {
                if (column[v]) {
                    row += '*';  // found a prime
                } else {
                    row += '-';  // found a composite
                }
            } else {
                row += ' ';  // not in sieveArray
            }
        }
        cerr << row << endl;
    }
}

template <>
void SieveTemplate<uint32_t>::printSieveArray() {
    // Print sieve array with same orientation as that in the complex plane.
    uint32_t columnMaxSize = 0;
    for (auto &column : sieveArray) {
        if (column.size() > columnMaxSize) {
            columnMaxSize = column.size();
        }
    }
    // Some type casting because subtraction confuses it.
    for (auto v = (int32_t)columnMaxSize - 1; v >= 0; v--) {
        string row;
        for (auto &column : sieveArray) {
            if (column.size() > v) {
                // printing entire block padded by 1
                // could use binary or hex or ints ...
                // bitset<32> b(sieveArray[u][v]);
                // row += b.to_string('-', '*');
                // row += ' ';
                stringstream stream;
                stream << setfill('0') << setw(8) << hex << column[v];
                string result(stream.str());
                row += result;
                row += ' ';
            } else {  // [u][v] index not in sieveArray
                row += string(9, ' ');
            }
        }
        cerr << row << endl;
    }
}

template <>
bool SieveTemplate<bool>::getSieveArrayValue(uint32_t u, uint32_t v) {
    return sieveArray.at(u).at(v);
}

template <>
uint32_t SieveTemplate<uint32_t>::getSieveArrayValue(uint32_t u, uint32_t v) {
    return sieveArray.at(u).at(v);
}




// Other useful general purpose functions.

// Integer square root.
uint32_t isqrt(uint64_t n) {
    uint64_t x, y;
    x = n;
    y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return (uint32_t)x;
}

// Positive remainder.
uint32_t mod(int64_t k, uint32_t m) {
    int32_t r = k % m;
    if (r < 0) {
        r += m;
    }
    return r;
}