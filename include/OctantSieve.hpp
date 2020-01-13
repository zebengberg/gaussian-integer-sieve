#pragma once
#include "BaseSieve.hpp"
using namespace std;

class OctantSieve : public SieveTemplate<bool> {
private:
    const uint64_t x;

public:
    // Using an initializer list and calling SieveTemplate constructor to set maxNorm
    explicit OctantSieve(uint64_t x, bool display = true) : SieveTemplate<bool>(x, display), x(x) {}

    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    uint64_t getCountBigPrimes() override;
};