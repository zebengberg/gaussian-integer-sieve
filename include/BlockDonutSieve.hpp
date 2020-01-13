#pragma once
#include "BaseSieve.hpp"
using namespace std;

class BlockDonutSieve : public SieveTemplate<uint32_t> {
private:
    uint32_t x, y, dx, dy;
    const uint32_t dStart[10];  // used to start d during crossOffMultiples()
    const uint32_t gapDonut[10][10];  // used to jump d during crossOffMultiples()
    const unsigned char bitDonut[10][10];  // used to compress a gint into a bit position
    const int32_t realPartDecompress[32];  // used to decompress a bit position into a gint
    const int32_t imagPartDecompress[32];

public:
    BlockDonutSieve(uint32_t, uint32_t, uint32_t, uint32_t, bool = true);
    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    uint64_t getCountBigPrimes() override;
};