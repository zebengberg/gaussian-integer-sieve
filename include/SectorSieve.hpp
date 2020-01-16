#pragma once
#include "BaseSieve.hpp"
#include <cmath>
using namespace std;

class SectorSieve : public SieveTemplate<bool> {
private:
    uint64_t x;
    double alpha, beta;
    const double tolerance = .00000000001;
    vector<uint32_t> heightShifts; // might not need

public:
    SectorSieve(uint64_t, double = 0, double = M_PI_4, bool = true);

    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    uint64_t getCountBigPrimes() override;
};