#pragma once
#include <vector>
using namespace std;

vector<pair<long, long>> gPrimes(long);
vector<pair<long, long>> gPrimesSegment(long, long, long);
unsigned long gPrimesCount(long);
unsigned long gPrimesSegmentCount(long, long, long);
pair<long *, unsigned long> gPrimesArray(long);  // holds pointer to array and size of array