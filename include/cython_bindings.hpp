#pragma once
#include <vector>
using namespace std;

vector<pair<uint32_t, uint32_t>> gPrimes(uint64_t);
vector<pair<uint32_t, uint32_t>> gPrimesSegment(uint32_t, uint32_t, uint32_t);
uint64_t gPrimesCount(uint64_t);
uint64_t gPrimesSegmentCount(uint32_t, uint32_t, uint32_t);
pair<uint32_t *, uint64_t> gPrimesAsArray(uint64_t);  // holds pointer to array of uint32_t and size of array