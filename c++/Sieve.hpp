#pragma once
#include <vector>
using namespace std;


// Gaussian integer struct.
struct gint {
    // Members are public by default in a struct.
    long a, b;
    gint(long a, long b) { this->a = a; this->b = b; }
    long norm() { return a * a + b * b; }
    gint flip() { return gint{b, a}; }  // clang likes curly brace list initialization
};


// Abstract base class to be used in various sieving implementations.
class SieveBase {
protected:
    long maxNorm;
    double progress;
    double totalProgress;
    vector<gint> smallPrimes;
    vector<gint> bigPrimes;

public:
    explicit SieveBase(long);  // constructor
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
// Cannot have an abstract base class template, so this derived class is a template.
template <typename T>
class SieveTemplate : public SieveBase {
protected:
    vector<vector<T>> sieveArray;  // T will be bool or unsigned int

public:
    explicit SieveTemplate(long maxNorm) : SieveBase(maxNorm) {};
    void printSieveArrayInfo();
    T getSieveArrayValue(long u, long v);
};


// Useful functions
long isqrt(long);
vector<gint> readPrimesFromFile(long);