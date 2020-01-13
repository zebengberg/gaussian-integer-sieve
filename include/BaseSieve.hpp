#pragma once
#include <vector>
using namespace std;



// Gaussian integer struct with 32-bit integer coordinates.
struct gint {
    // Members are public by default in a struct.
    uint32_t a, b;
    gint(uint32_t a, uint32_t b) { this->a = a; this->b = b; }
    uint64_t norm() { return a * a + b * b; }
    gint flip() { return gint{b, a}; }  // clang likes curly brace list initialization
    pair<uint32_t, uint32_t> asPair() {return pair<uint32_t, uint32_t> {a, b}; }
    friend bool operator < (gint g1, gint g2) { return g1.norm() < g2.norm(); }
};


// Abstract base class to be used in various sieving implementations.
class SieveBase {
protected:
    const uint64_t maxNorm;
    double progress;
    double totalProgress;
    bool verbose;
    vector<gint> smallPrimes;
    vector<gint> bigPrimes;

public:
    explicit SieveBase(uint64_t, bool);  // constructor; will be called in derived classes
    void setSmallPrimesFromFile();
    void sieve();  // crossing off all multiples of small primes
    void printProgress(gint);
    void sortBigPrimes();
    void printBigPrimes();
    void writeBigPrimesToFile();
    void run();  // run necessary sieve methods; does not gather big primes
    vector<gint> getBigPrimes();  // return the big primes after run()

    // Virtual methods to be implemented in derived classes
    virtual void setSmallPrimes() = 0;
    virtual void setSieveArray() = 0;
    virtual void crossOffMultiples(gint) = 0;
    virtual void setBigPrimes() = 0;  // results of sieve
    virtual uint64_t getCountBigPrimes() = 0;
};


// Derived from SieveBase and parameterized with sieveArray entry type.
// Cannot have an abstract base class template, so this derived class is a template.
template <typename T>
class SieveTemplate : public SieveBase {
protected:
    vector<vector<T>> sieveArray;  // T will be bool or unsigned int

public:
    explicit SieveTemplate(uint64_t maxNorm, bool display) : SieveBase(maxNorm, display) {};
    void printSieveArrayInfo();
    void printSieveArray();
    T getSieveArrayValue(uint32_t u, uint32_t v);
};


// Useful library-style functions
uint32_t isqrt(uint64_t);