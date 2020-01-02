#pragma once
#include "Sieve.hpp"
using namespace std;

class SegmentedSieve : public SieveTemplate<bool> {
private:
    long x, y, z;

public:
    SegmentedSieve(long, long, long, bool = true);
    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
};