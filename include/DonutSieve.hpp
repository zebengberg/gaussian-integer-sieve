#pragma once
#include "BaseSieve.hpp"
#include "SegmentedSieve.hpp"
using namespace std;

class DonutSieve : public SieveTemplate<unsigned int> {
private:
    long x;
    int dStart[10];  // used to start d during crossOffMultiples()
    int gapDonut[10][10];  // used to jump d during crossOffMultiples()
    unsigned int bitDonut[10][10];  // used to compress a gint into a bit position
    int realPartDecompress[32];  // used to decompress a bit position into a gint
    int imagPartDecompress[32];

public:
    explicit DonutSieve(long, bool = true);  // default values must be set in header
    void setFalse(long, long);
    void setTrue(long, long);
    static void printDonutArrays();  // only needed to help write the source cpp file
    // overriding virtual methods
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    unsigned long getCountBigPrimes() override;

};