#include <iostream>
#include <cmath>
#include "gint.hpp"
using namespace std;

class SegmentedSieve {
private:
    long x, y, z;
    vector<vector<bool>> sieveArray;
    vector<Gint> primes;  // to hold primes after sieving

public:
    SegmentedSieve(long, long, long)
};

SegmentedSieve::SegmentedSieve(long x, long y, long z) {
    this->x = x;  // x-coordinate of lower left-hand corner of segment block
    this->y = y;  // y-coordinate of lower left-hand corner of segment block
    this->z = z;  // side length of segment block
}