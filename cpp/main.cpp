#include "main.hpp"
#include "QuadrantSieve.hpp"
#include "DonutSieve.hpp"
#include <vector>
#include <cmath>
//using namespace std;


//int main(int argc, const char* argv[]){
//    if (argc != 2) {
//        cerr << "Specify the sieve array norm bound." << endl;
//        return 1;
//    }
//    long x = strtol(argv[1], nullptr, 10);  // convert command line input to long


//int main() {
//    long x = pow(10, 9);
//    long y = pow(10, 8);
//    long z = 100;
//    SegmentedSieve s(x, y, z);
//    s.run();
//    s.countBigPrimes();
//    return 0;
//}

int main() {
    long x = 10000;
    QuadrantSieve s(x);
    s.run();
    s.sortBigPrimes();
    s.printBigPrimes();
    // s.printSieveArray();
    // s.countBigPrimes();
    return 0;
}


vector<pair<long, long>> gPrimes(long x) {
    // Show display if passed argument is large.
    bool display = x >= (long)pow(10, 7);
    DonutSieve s(x, display);
    vector<gint> gintP = s.run();
    vector<pair<long, long>> pairP;
    pairP.resize(gintP.size());
    transform(gintP.begin(), gintP.end(), pairP.begin(), [](gint g) { return g.asPair(); });  // lambda
    return pairP;
}


long gPrimesCount(long x) {
    // Show display if passed argument is large.
    bool display = x >= (long)pow(10, 7);
    DonutSieve s(x, display);
    vector<gint> gintP = s.run();
    return gintP.size();
}