#pragma once
#include <vector>
using namespace std;


// Gaussian integer struct.
struct gint {
    // Members are public by default in a struct.
    long a, b;
    gint(long a, long b) { this->a = a; this->b = b; }
    long norm() { return a * a + b * b; }
};


// Abstract base class to be used in various sieving implementations.
class SieveBase {
protected:
    long maxNorm;
    double progress = 0.0;
    double totalProgress;
    vector<gint> smallPrimes;
    vector<gint> bigPrimes;

public:
    virtual void setMemberVariables() = 0;
    virtual void setSmallPrimes() = 0;
    virtual void setSieveArray() = 0;
    virtual void crossOffMultiples(gint) = 0;
    void sieve();  // crossing off all multiples of small primes
    void printProgress(gint);
    virtual void setBigPrimes() = 0;  // results of sieve
    void printBigPrimes();
    void countBigPrimes();
    void writeBigPrimesToFile();
};

// Derived from SieveBase and parameterized with sieveArray entry type.
// Cannot have an abstract base class template.
template <typename T>
class SieveTemplate : public SieveBase {
protected:
    vector<vector<T>> sieveArray;  // T will be bool or int

public:
    void printSieveArrayInfo();
};


// Useful functions
long isqrt(long);
vector<gint> readPrimesFromFile(long);