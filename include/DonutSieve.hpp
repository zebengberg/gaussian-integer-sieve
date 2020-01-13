#pragma once
#include "BaseSieve.hpp"
#include "SegmentedSieve.hpp"
using namespace std;

class DonutSieve : public SieveTemplate<unsigned int> {
private:
    const uint64_t x;
    const uint32_t dStart[10];  // used to start d during crossOffMultiples()
    const uint32_t gapDonut[10][10];  // used to jump d during crossOffMultiples()
    const unsigned char bitDonut[10][10];  // used to compress a gint into a bit position
    const int32_t realPartDecompress[32];  // used to decompress a bit position into a gint
    const int32_t imagPartDecompress[32];

public:
    explicit DonutSieve(uint64_t, bool = true);  // default values must be set in header
    void setFalse(uint32_t, uint32_t);
    void setTrue(uint32_t, uint32_t);
    static void printDonutArrays();  // only needed to help write the source cpp file
    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    uint64_t getCountBigPrimes() override;

};