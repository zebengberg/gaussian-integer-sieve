#pragma once
#include "BaseSieve.hpp"
using namespace std;

class SegmentedSieve : public SieveTemplate<bool> {
private:
    uint32_t x, y, z, w;

public:
    SegmentedSieve(uint32_t, uint32_t, uint32_t, uint32_t, bool = true);
    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    uint64_t getCountBigPrimes() override;
};