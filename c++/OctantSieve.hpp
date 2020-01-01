#pragma once
#include "Sieve.hpp"
using namespace std;

class OctantSieve : public SieveTemplate<bool> {
private:
    long x;

public:
    // Calling SieveTemplate constructor to set maxNorm
    explicit OctantSieve(long x) : SieveTemplate<bool>(x) { this->x = x;}
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
};