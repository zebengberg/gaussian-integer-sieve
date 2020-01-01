#pragma once
#include "Sieve.hpp"


class QuadrantSieve : public SieveTemplate<bool> {
private:
    long x;

public:
    // Using an initializer list and calling SieveTemplate constructor to set maxNorm
    explicit QuadrantSieve(long x) : x(x), SieveTemplate<bool>(x) {}
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void crossOffMultiplesAlt(gint);  // old method; keeping for reference
    void setBigPrimes() override;
};