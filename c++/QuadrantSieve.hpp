#pragma once
#include "Sieve.hpp"


class QuadrantSieve : public SieveTemplate<bool> {
private:
    long x;

public:
    explicit QuadrantSieve(long x) { this->x = x;}
    void setMemberVariables() override;
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void crossOffMultiplesAlt(gint);  // old method; keeping for reference
    void setBigPrimes() override;
};