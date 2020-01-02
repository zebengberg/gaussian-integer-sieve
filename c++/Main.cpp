#include "QuadrantSieve.hpp"
//#include "OctantSieve.hpp"
//#include "SegmentedSieve.hpp"
//#include "DonutSieve.hpp"
#include <cmath>
using namespace std;


//int main(int argc, const char* argv[]){
//    if (argc != 2) {
//        cerr << "Specify the sieve array norm bound." << endl;
//        return 1;
//    }
//    long x = strtol(argv[1], nullptr, 10);  // convert command line input to long


//int main() {
//    long x = 10000000;
//    long y = 5000000;
//    long z = 1000;
//    SegmentedSieve s(x, y, z);
//    s.setSieveArray();
//    s.printSieveArrayInfo();
//    s.setSmallPrimesSegmented();
//    s.sieve();
//    s.setBigPrimes();
//    s.countBigPrimes();
//    s.printSieveArray();
//    return 0;
//}

int main() {
    long x = 10000000;
    QuadrantSieve s(x, false);
    s.setSmallPrimes();
    s.setSieveArray();
    s.sieve();
    //s.printSieveArray();
    s.setBigPrimes();
    s.countBigPrimes();
    return 0;
}

// TODO: run code inspect