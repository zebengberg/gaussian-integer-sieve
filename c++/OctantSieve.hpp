#pragma once
#include "Sieve.hpp"
using namespace std;

class OctantSieve : public SieveTemplate<bool> {
private:
    long x;

public:
    // Using an initializer list and calling SieveTemplate constructor to set maxNorm
    explicit OctantSieve(long x, bool display = true) : x(x), SieveTemplate<bool>(x, display) {}

    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
};