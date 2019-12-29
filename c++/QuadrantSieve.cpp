#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include "gint.hpp"
using namespace std;

// Read small primes to be used for sieving from existing file.
// Getting all primes from file small_primes.txt with norm up to x.
vector<Gint> readPrimesFromFile(long x) {
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
    vector<Gint> smallPrimes;  // the list of small primes to be populated
    while (a * a + b * b <= x) {
        Gint p(a, b);  // creating prime struct
        smallPrimes.push_back(p);  // putting it into vector P

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

class QuadrantSieve{
private:
    long x;
    double progress;
    double totalProgress;
    vector<vector<bool>> sieveArray;
    vector<Gint> primes;  // to hold primes after sieving

public:
    explicit QuadrantSieve(long);  // constructor
    void initializeSieveArray();
    void crossOffMultiplesSkew(Gint);
    void crossOffMultiplesRect(Gint);
    void getPrimes();
    void printPrimes();
    void writePrimesToFile();
    void printProgress(Gint);  // update after each sieved prime
};

QuadrantSieve::QuadrantSieve(long x) {
    this->x = x;
    progress = 0.0;
    totalProgress = log(log(x)) - log(2.0);
}

void QuadrantSieve::initializeSieveArray() {
    // Lots of controversy about vector<bool>; google it.
    // Each boolean value is stored as a single bit but pointers are not available.
    // sieveArray holds values for Gint's with a, b >= 0 and a^2 + b^2 <= x.
    cout << "Building sieve array..." << endl;
    for (long a = 0; a <= isqrt(x); a++) {
        long b = isqrt(x - a * a) + 1;
        vector<bool> column((unsigned long)b, true);  // Create a vector of size b with all values true.
        sieveArray.push_back(column);
    }
    sieveArray[0][0] = false;  // 0 is not prime
    sieveArray[1][0] = false;  // 1 is not prime
    sieveArray[0][1] = false;  // i is not prime

    double size = 3.1415 * double(x) / 4;  // area of the quarter circle
    size /= 8;  // 8 bits per byte
    size /= 1000000000;  // convert to GB
    cout << "Sieve array approximate memory use: " << size  << "GB" << endl;
}

void QuadrantSieve::crossOffMultiplesSkew(Gint g) {
    for (long c = 1; c <= isqrt(x); c++) {
        long u = c * g.a();  // will become negative
        long v = c * g.b();
        for (long d = 0; d <= isqrt(x / g.norm() - c * c); d++) {
            if (u > 0) {  // inexpensive; could rewrite to avoid this, but would be less understandable
                sieveArray[u][v] = false;
            } else {  // multiply u + vi by -i to bring it back to the first quadrant
                sieveArray[v][-u] = false;
            }
            u -= g.b();
            v += g.a();
        }
    }
    sieveArray[g.a()][g.b()] = true;  // crossed this off; need to remark it as prime
}

void QuadrantSieve::crossOffMultiplesRect(Gint g) {
    long p, subgroupSize;
    if (g.b()) {  // degree 1 primes
        p = g.norm();
        subgroupSize = p;
    } else {  // degree 2 primes
        p = g.a();
        subgroupSize = 1;  // trivial subgroup
    }
    for (long i = 0; i < subgroupSize; i++) {
        // s + it is an element of the additive subgroup < a + bi > in Z[i] / pZ[i]
        long s = i * g.a() % p;
        long t = i * g.b() % p;
        for (long u = s; u <= isqrt(x); u += p) {
            for (long v = t; v <= isqrt(x - u * u); v += p) {
                sieveArray[u][v] = false;
            }
        }
    }
    sieveArray[g.a()][g.b()] = true;  // crossed this off; need to remark it as prime
}

void QuadrantSieve::getPrimes() {
    cout << "Gathering primes after sieve..." << endl;
    for (long a = 1; a <= isqrt(x); a++) {  // Want to stay off imaginary line
        for (long b = 0; b <= isqrt(x - a * a); b++) {
            if (sieveArray[a][b]) {
                Gint p(a, b);
                primes.push_back(p);
            }
        }
    }
}

void QuadrantSieve::printPrimes() {
    for (Gint p : primes) {
        cout << p.a() << "  " << p.b() << "  " << p.norm() << endl;
    }
    cout << "Total number of primes: " << primes.size() << endl;
}

void QuadrantSieve::writePrimesToFile() {
    ofstream f;
    f.open("../data/cpp_primes.csv");
    for (Gint p : primes) {
        f << p.a() << "," << p.b() << "," << p.norm() << endl;
    }
    f.close();
}

void QuadrantSieve::printProgress(Gint g) {
    int barSize = 70;
    progress += 1.0 / double(g.norm());
    int barPos = int(double(barSize) * progress / totalProgress);
    cout << "[";
    for (int i = 0; i < barPos; i++) { cout << "|"; }
    for (int i = barPos; i <= barSize; i++) { cout << " "; }
    cout << "]\r";
    cout.flush();
}


int main(int argc, const char* argv[]){
    if (argc != 2) {
        cerr << "Specify the sieve array norm bound." << endl;
        return 1;
    }
    long x = strtol(argv[1], nullptr, 10);  // convert command line input to long
    QuadrantSieve qs(x);
    qs.initializeSieveArray();
    vector<Gint> smallPrimes = readPrimesFromFile(isqrt(x));
    cout << "Starting to sieve..." << endl;
    auto startTime = chrono::high_resolution_clock::now();
    for (Gint p : smallPrimes) {
        qs.crossOffMultiplesRect(p);
        qs.printProgress(p);
    }
    auto endTime = chrono::high_resolution_clock::now();
    auto totalTime = chrono::duration_cast<chrono::seconds>(endTime - startTime);
    cout << "Total time for sieving: " << totalTime.count() << " seconds" << endl;
    qs.getPrimes();
    qs.printPrimes();
    //sieveArray.writePrimesToFile();
    return 0;
}
