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
class Sieve {
protected:
    long maxNorm;
    double progress = 0.0;
    double totalProgress;
    vector<vector<bool>> sieveArray;
    vector<gint> smallPrimes;
    vector<gint> bigPrimes;

public:
    virtual void setMemberVariables() = 0;
    virtual void setSmallPrimes() = 0;
    virtual void setSieveArray() = 0;
    void printSieveArrayInfo();
    virtual void crossOffMultiples(gint) = 0;
    void sieve();
    void printProgress(gint);
    virtual void setBigPrimes() = 0;  // results of sieve
    void printBigPrimes();
    void countBigPrimes();
    void writeBigPrimesToFile();
};

// Useful functions
long isqrt(long);
vector<gint> readPrimesFromFile(long);