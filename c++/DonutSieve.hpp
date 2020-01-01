#pragma once
#include "Sieve.hpp"
using namespace std;

class DonutSieve : public SieveTemplate<unsigned int> {
private:
    long x;
    int dStart[10];  // used to start d during crossOffMultiples()
    int gapDonut[10][10];  // used to jump d during crossOffMultiples()
    int bitDonut[10][10];  // used to compress a gint into a bit position
    int realPartDecompress[32];  // used to decompress a bit position into a gint
    int imagPartDecompress[32];

public:
    explicit DonutSieve(long);
    void setSmallPrimes() override;
    void setSieveArray() override;
    void setFalse(long, long);
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    void printDonut();
};