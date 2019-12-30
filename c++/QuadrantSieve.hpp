#pragma once
#include "Sieve.hpp"


class QuadrantSieve : public Sieve {
private:
    long x;

public:
    explicit QuadrantSieve(long x) { this->x = x;}
    void setMemberVariables() override;
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
};