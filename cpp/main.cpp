#include "main.hpp"
#include "QuadrantSieve.hpp"
#include "OctantSieve.hpp"
#include "DonutSieve.hpp"
#include "SegmentedSieve.hpp"
#include <cmath>
//using namespace std;


//int main(int argc, const char* argv[]){
//    if (argc != 2) {
//        cerr << "Specify the sieve array norm bound." << endl;
//        return 1;
//    }
//    long x = strtol(argv[1], nullptr, 10);  // convert command line input to long



int main() {
    long x = 20000;
    DonutSieve s(x);
    s.run();
    s.countBigPrimes();
    s.printBigPrimes();
    // s.printSieveArray();
    return 0;
}


vector<pair<long, long>> gPrimes(long x) {
    // Show display if passed argument is large.
    bool display = x >= (long)pow(10, 7);
    DonutSieve s(x, display);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    vector<pair<long, long>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    return pairP;
}

unsigned long gPrimesCount(long x) {
    // Show display if passed argument is large.
    bool display = x >= (long)pow(10, 7);
    DonutSieve s(x, display);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    return 4 * gintP.size();  // four quadrants
}

vector<pair<long, long>> gPrimesSegment(long x, long y, long z) {
    // Show display if passed argument is large.
    bool display = z >= (long)pow(10, 5);
    SegmentedSieve s(x, y, z, display);
    s.run();
    vector<gint> gintP = s.getBigPrimes();
    vector<pair<long, long>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    return pairP;
}