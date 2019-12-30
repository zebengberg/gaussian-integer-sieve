#pragma once
#include "Sieve.hpp"
using namespace std;

class SegmentedSieve : public Sieve {
private:
    long x, y, z;


public:
    SegmentedSieve(long, long, long);
    void setMemberVariables() override;
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;

    void printSieveArray();
};