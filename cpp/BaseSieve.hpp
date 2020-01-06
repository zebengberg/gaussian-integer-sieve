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
    pair<long, long> asPair() {return pair<long, long> {a, b}; }
    friend bool operator < (gint g1, gint g2) { return g1.norm() < g2.norm(); }
};


// Abstract base class to be used in various sieving implementations.
class SieveBase {
protected:
    long maxNorm;
    double progress;
    double totalProgress;
    bool display;
    vector<gint> smallPrimes;
    vector<gint> bigPrimes;

public:
    explicit SieveBase(long, bool);  // constructor; will be called in derived classes
    void setSmallPrimesFromFile();
    void sieve();  // crossing off all multiples of small primes
    void printProgress(gint);
    void sortBigPrimes();
    void printBigPrimes();
    void countBigPrimes();
    void writeBigPrimesToFile();
    vector<gint> run();  // run everything and return sorted bigPrimes

    // Virtual methods to be implemented in derived classes
    virtual void setSmallPrimes() = 0;
    virtual void setSieveArray() = 0;
    virtual void crossOffMultiples(gint) = 0;
    virtual void setBigPrimes() = 0;  // results of sieve
};


// Derived from SieveBase and parameterized with sieveArray entry type.
// Cannot have an abstract base class template, so this derived class is a template.
template <typename T>
class SieveTemplate : public SieveBase {
protected:
    vector<vector<T>> sieveArray;  // T will be bool or unsigned int

public:
    explicit SieveTemplate(long maxNorm, bool display) : SieveBase(maxNorm, display) {};
    void printSieveArrayInfo();
    void printSieveArray();
    T getSieveArrayValue(unsigned long u, unsigned long v);
};


// Useful functions
long isqrt(long);
vector<gint> readPrimesFromFile(long);