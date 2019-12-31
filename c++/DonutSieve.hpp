#pragma once
#include "Sieve.hpp"
using namespace std;

class DonutSieve : public Sieve {
private:
    long x;

public:
    explicit DonutSieve(long x) { this->x = x; }
    void setMemberVariables() override;
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    gint extractGint(int, int);
    void printDonut();
};