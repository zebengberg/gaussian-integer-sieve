#pragma once
#include "Sieve.hpp"
using namespace std;



class OctantSieve : public Sieve {
private:
    long x;

public:
    explicit OctantSieve(long x) { this-> x = x; }
    void setMemberVariables() override;
    void setSmallPrimes() override;
    void setSieveArray() override;
    void crossOffMultiples(gint) override;
    void setBigPrimes() override;
};