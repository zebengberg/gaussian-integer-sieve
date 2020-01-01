#pragma once
#include "Sieve.hpp"
using namespace std;

class SegmentedSieve : public SieveTemplate<bool> {
private:
    long x, y, z;

public:
    SegmentedSieve(long, long, long);
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
};