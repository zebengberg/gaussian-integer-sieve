#pragma once
#include <vector>
using namespace std;

// Returning containers which get slowly copied into python structures.
vector<pair<uint32_t, uint32_t>> gPrimesToNorm(uint64_t);
vector<pair<uint32_t, uint32_t>> gPrimesInSector(uint64_t, double, double);
vector<pair<uint32_t, uint32_t>> gPrimesInBlock(uint32_t, uint32_t, uint32_t, uint32_t);

// Returning counts only.
uint64_t gPrimesToNormCount(uint64_t);
uint64_t gPrimesInSectorCount(uint64_t, double, double);
uint64_t gPrimesInBlockCount(uint32_t, uint32_t, uint32_t, uint32_t);

// Returning pointers that can be shared with numpy to build np.arrays.
// The first element of the pair holds a pointer to an array of uint32_t ints,
// and the second element of the pair holds the size of the array.
// Real parts and imaginary parts are concatenated together into the array;
// they can be untangled in numpy as needed.
pair<uint32_t *, uint64_t> gPrimesToNormAsArray(uint64_t);
pair<uint32_t *, uint64_t> gPrimesInSectorAsArray(uint64_t, double, double);
pair<uint32_t *, uint64_t> gPrimesInBlockAsArray(uint32_t, uint32_t, uint32_t, uint32_t);

// Histogram of angles of primes to norm.
vector<uint64_t> angularDistribution(uint64_t, uint32_t);
int32_t * sectorRace(uint64_t, double, double, double, double, uint32_t);