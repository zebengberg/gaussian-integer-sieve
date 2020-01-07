#pragma once
#include "BaseSieve.hpp"


class QuadrantSieve : public SieveTemplate<bool> {
private:
    long x;

public:
    // Using an initializer list and calling SieveTemplate constructor to set maxNorm
    explicit QuadrantSieve(long x, bool display = true) : SieveTemplate<bool>(x, display), x(x) {}

    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void crossOffMultiplesAlt(gint);  // old method; keeping for reference
    void setBigPrimes() override;
    unsigned long getCountBigPrimes() override;
};