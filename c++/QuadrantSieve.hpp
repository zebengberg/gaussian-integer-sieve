#pragma once
#include "Sieve.hpp"


class QuadrantSieve : public SieveTemplate<bool> {
private:
    long x;

public:
    // Calling SieveTemplate constructor to set maxNorm
    explicit QuadrantSieve(long x) : SieveTemplate<bool>(x) { this->x = x;}
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void crossOffMultiplesAlt(gint);  // old method; keeping for reference
    void setBigPrimes() override;
};