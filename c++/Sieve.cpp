#include "Sieve.hpp"
#include <iostream>
#include <fstream>
#include <cmath>



template <typename T>
void Sieve<T>::printSieveArrayInfo() {
    unsigned long totalSize = sizeof(sieveArray);
    for (const auto& column : sieveArray) {
        totalSize += sizeof(column);
        if (is_same<T, bool>::value) {
            totalSize += column.capacity() / 8;  // each bool stored as a bit
        } else {
            totalSize += column.capacity() * sizeof(T);
        }
    }
    totalSize /= pow(10, 6);  // convert to MB
    cout << "Sieve array approximate memory use: " << totalSize  << "MB" << endl;
}

template <typename T>
void Sieve<T>::sieve() {
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

template <typename T>
void Sieve<T>::printProgress(gint g) {
    int barSize = 80;
    progress += 1.0 / double(g.norm());
    int barPos = int(double(barSize) * progress / totalProgress);
    cout << "[";
    for (int i = 0; i < barPos; i++) { cout << "|"; }
    for (int i = barPos; i <= barSize; i++) { cout << " "; }
    cout << "]\r";
    cout.flush();
}

template <typename T>
void Sieve<T>::printBigPrimes() {
    for (gint g : bigPrimes) {
        cout << g.a << "  " << g.b << "  " << g.norm() << endl;
    }
    Sieve::countBigPrimes();
}

template <typename T>
void Sieve<T>::countBigPrimes() {
    cout << "Total number of primes, including associates: " << 4 * bigPrimes.size() << endl;
}

template <typename T>
void Sieve<T>::writeBigPrimesToFile() {
    ofstream f;
    f.open("../data/cpp_primes.csv");
    for (gint g : bigPrimes) {
        f << g.a << ' ' << g.b << endl;
    }
    f.close();
}

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


// Read small primes to be used for sieving from existing file.
// Getting all primes from file small_primes.txt with norm up to x.
vector<gint> readPrimesFromFile(long x) {
    // If x is too small, then we are not actually taking any primes.
    if (x < 2) {
        cerr << "No primes read and no sieving needed.";
        exit(1);
    }
    // small_primes.txt contains a list of Gaussian primes sorted by norm.
    ifstream f;
    f.open("../data/small_primes.txt");
    if (!f) {
        cerr << "Unable to open file small_primes.txt";
        exit(1);   // call system to stop
    } else {
        cout << "Reading small primes from file." << endl;
    }

    long a, b;
    f >> a >> b;  // token-based parsing
    vector<gint> smallPrimes;  // the list of small primes to be populated
    while (a * a + b * b <= x) {
        gint g(a, b);  // creating prime struct
        smallPrimes.push_back(g);  // putting it into vector P

        f >> a >> b;  // reading next pair a, b from file f
        // If we've gotten to the end of the file, then we don't have enough precomputed primes.
        if (f.eof()) {
            cerr << "Not enough primes in small_primes.txt";
            f.close();
            exit(1);
        }
    }
    f.close();
    return smallPrimes;
}
