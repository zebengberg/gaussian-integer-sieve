#pragma once
#include "Sieve.hpp"
using namespace std;

class OctantSieve : public SieveTemplate<bool> {
private:
    long x;

public:
    // Using an initializer list and calling SieveTemplate constructor to set maxNorm
    explicit OctantSieve(long x) : x(x), SieveTemplate<bool>(x) {}
    // overriding virtual methods
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
};