#include "Sieve.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

// Will call this constructor from derived classes.
SieveBase::SieveBase(long maxNorm) : maxNorm(maxNorm) {  // initializer list
    progress = 0.0;
    // Using PNT to make approximate progress bar.
    totalProgress = log(log(maxNorm)) - log(2.0);
}

void SieveBase::sieve() {
    cout << "Starting to sieve..." << endl;
    auto startTime = chrono::high_resolution_clock::now();
    for (gint g : smallPrimes) {
        crossOffMultiples(g);
        printProgress(g);
    }
    auto endTime = chrono::high_resolution_clock::now();
    auto totalTime = chrono::duration_cast<chrono::seconds>(endTime - startTime);
    cout << "Total time for sieving: " << totalTime.count() << " seconds" << endl;
}

void SieveBase::printProgress(gint g) {
    int barSize = 80;
    progress += 1.0 / double(g.norm());
    int barPos = int(double(barSize) * progress / totalProgress);
    cout << "[";
    for (int i = 0; i < barPos; i++) { cout << "|"; }
    for (int i = barPos; i <= barSize; i++) { cout << " "; }
    cout << "]\r";
    cout.flush();
}

vector<gint> SieveBase::getBigPrimes() { return bigPrimes; }

void SieveBase::printBigPrimes() {
    for (gint g : bigPrimes) {
        cout << g.a << "  " << g.b << "  " << g.norm() << endl;
    }
    SieveBase::countBigPrimes();
}

void SieveBase::countBigPrimes() {
    cout << "Total number of primes, including associates: " << 4 * bigPrimes.size() << endl;
}

void SieveBase::writeBigPrimesToFile() {
    ofstream f;
    f.open("../data/cpp_primes.csv");
    for (gint g : bigPrimes) {
        f << g.a << ' ' << g.b << endl;
    }
    f.close();
}

// Getting all primes from file small_primes.txt with norm up to maxNorm variable.
void SieveBase::setSmallPrimes() {
    // If maxNorm is too small, then we are not actually taking any primes.
    if (maxNorm < 2) {
        cerr << "No primes read and no sieving needed.";
        exit(1);
    }
    // small_primes.txt contains a list of Gaussian primes sorted by norm.
    ifstream f;
    f.open("../data/small_primes.txt");
    if (!f) {
        cerr << "Unable to open file small_primes.txt";
        exit(1);
    } else {
        cout << "Reading small primes from file." << endl;
    }

    long a, b;
    f >> a >> b;  // token-based parsing
    gint g(a, b);
    // Need primes up to square root of the maxNorm
    while (g.norm() <= isqrt(maxNorm)) {
        smallPrimes.push_back(g);
        f >> a >> b;  // reading next pair a, b from file f
        g = gint(a, b);
        // If we get to the end of the file, then we don't have enough precomputed primes.
        if (f.eof()) {
            cerr << "Not enough primes in small_primes.txt!" << endl;
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
    unsigned long totalSize = sizeof(sieveArray);
    for (const auto& column : sieveArray) {
        totalSize += sizeof(column) + column.capacity() / 8;  // each bool stored as a bit
    }
    totalSize /= pow(10, 6);  // convert to MB
    cout << "Sieve array approximate memory use: " << totalSize  << "MB" << endl;
}

template <>
void SieveTemplate<unsigned int>::printSieveArrayInfo() {
    unsigned long totalSize = sizeof(sieveArray);
    for (const auto& column : sieveArray) {
        totalSize += sizeof(column) + column.capacity() * sizeof(unsigned int);
    }
    totalSize /= pow(10, 6);  // convert to MB
    cout << "Sieve array approximate memory use: " << totalSize  << "MB" << endl;
}

template <>
void SieveTemplate<bool>::printSieveArray() {
    // Print sieve array with same orientation as that in the complex plane.
    long columnMaxSize = 0;
    for (const auto& column : sieveArray) {
        if (column.size() > columnMaxSize) {
            columnMaxSize = column.size();
        }
    }
    for (long v = columnMaxSize - 1; v >=0; v--) {
        string row;
        for (long u = 0; u < sieveArray.size(); u++) {
            if (sieveArray[u].size() > v) {
                if (sieveArray[u][v]) {
                    row += '*';  // found a prime
                } else {
                    row += '-';  // found a composite
                }
            } else {
                row += ' ';  // not in sieveArray
            }
        }
        cout << row << endl;
    }
}

template <>
void SieveTemplate<unsigned int>::printSieveArray() {
    // Print sieve array with same orientation as that in the complex plane.
    long columnMaxSize = 0;
    for (const auto& column : sieveArray) {
        if (column.size() > columnMaxSize) {
            columnMaxSize = column.size();
        }
    }
    for (long v = columnMaxSize - 1; v >=0; v--) {
        string row;
        for (long u = 0; u < sieveArray.size(); u++) {
            if (sieveArray[u].size() > v) {
                // printing entire block padded by 1
                // could use binary or hex or ints ...
                // bitset<32> b(sieveArray[u][v]);
                // row += b.to_string('-', '*');
                // row += ' ';
                stringstream stream;
                stream << setfill('0') << setw(8) << hex << sieveArray[u][v];
                string result(stream.str());
                row += result;
                row += ' ';
            } else {  // [u][v] index not in sieveArray
                row += string(9, ' ');
            }
        }
        cout << row << endl;
    }
}

template <>
bool SieveTemplate<bool>::getSieveArrayValue(long u, long v) { return sieveArray.at(u).at(v); }

template <>
unsigned int SieveTemplate<unsigned int>::getSieveArrayValue(long u, long v) { return sieveArray.at(u).at(v); }




// Other useful general purpose functions.

// Integer square root.
long isqrt(long n) {
    long x, y;
    x = n;
    y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

