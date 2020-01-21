#pragma once
#include "BaseSieve.hpp"
#include <cmath>
using namespace std;

class SectorSieve : public SieveTemplate<bool> {
private:
    uint64_t x;
    double alpha, beta;
    // Only need tolerance when alpha or beta is close to rational multiple of pi.
    const double tolerance = pow(10, -10);
    vector<int32_t> heightShifts;

public:
    SectorSieve(uint64_t, double, double, bool = true);

    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    uint64_t getCountBigPrimes() override;
};