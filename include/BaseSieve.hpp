#pragma once
#include <vector>
#include <cmath>
using namespace std;



// Gaussian integer struct with 32-bit integer coordinates.
struct gint {
    // Members are public by default in a struct.
    int32_t a, b;
    gint(int32_t a, int32_t b) { this->a = a; this->b = b; }
    uint64_t norm() { return uint64_t(a) * uint64_t(a) + uint64_t(b) * uint64_t(b); }
    double arg() { return atan2(b, a); }
    gint flip() { return gint{b, a}; }  // clang likes curly brace list initialization
    pair<int32_t, int32_t> asPair() {return pair<int32_t, int32_t> {a, b}; }
    // Gives a total ordering; every two distinct gints are related by <.
    friend bool operator < (gint g1, gint g2) {
        uint64_t n1 = g1.norm();
        uint64_t n2 = g2.norm();
        if (n1 == n2) {
            return g1.a > g2.a;  // reverse lexicographic; larger real part is smaller
        } else {
            return n1 < n2;
        }
    }
    friend bool operator == (gint g1, gint g2) {
        return (g1.a == g2.a) && (g1.b == g2.b);
    }
    friend gint operator + (gint g1, gint g2) {
        return gint(g1.a + g2.a, g1.b + g2.b);
    }
};


// Abstract base class to be used in various sieving implementations.
class SieveBase {
protected:
    const uint64_t maxNorm;
    double progress;
    double totalProgress;
    uint32_t discreteProgress;
    bool verbose;
    vector<gint> smallPrimes;
    vector<gint> bigPrimes;

public:
    explicit SieveBase(uint64_t, bool);  // constructor; will be called in derived classes
    void setSmallPrimesFromFile();
    void setSmallPrimesFromReference(const vector<gint>&);
    void sieve();  // crossing off all multiples of small primes
    void printProgress(gint);
    void sortBigPrimes();
    void printBigPrimes();
    void writeBigPrimesToFile();
    void run();  // run necessary sieve methods; does not gather big primes
    vector<gint> getBigPrimes(bool = true);  // return the big primes after run()

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
    vector<vector<T>> getSieveArray();
};


// Useful library-style functions
uint32_t isqrt(uint64_t);
uint32_t mod(int64_t, uint32_t);