#pragma once
#include "BaseSieve.hpp"


class QuadrantSieve : public SieveTemplate<bool> {
private:
    const uint64_t x;

public:
    // Using an initializer list and calling SieveTemplate constructor to set maxNorm
    explicit QuadrantSieve(uint64_t x, bool verbose = true) : SieveTemplate<bool>(x, verbose), x(x) {}

    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    uint64_t getCountBigPrimes() override;
};