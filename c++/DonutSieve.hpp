#pragma once
#include "Sieve.hpp"
using namespace std;

class DonutSieve : public SieveTemplate<unsigned int> {
private:
    long x;
    int dStart[10];  // used to start d during crossOffMultiples()
    int donut[10][10];  // used to jump d during crossOffMultiples()

public:
    explicit DonutSieve(long);
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
    //gint extractGint(int, int);
    void printDonut();
};