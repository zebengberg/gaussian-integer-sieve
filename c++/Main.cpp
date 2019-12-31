#include "QuadrantSieve.hpp"
//#include "OctantSieve.hpp"
//#include "SegmentedSieve.hpp"
//#include "DonutSieve.hpp"
#include <iostream>
using namespace std;



//int main(int argc, const char* argv[]){
//    if (argc != 2) {
//        cerr << "Specify the sieve array norm bound." << endl;
//        return 1;
//    }
//    long x = strtol(argv[1], nullptr, 10);  // convert command line input to long


//int main() {
//    long x = 30000;
//    long y = 10000;
//    long z = 150;
//    SegmentedSieve s(x, y, z);
//    s.setMemberVariables();
//    s.setSieveArray();
//    s.printSieveArrayInfo();
//    s.setSmallPrimes();
//    s.sieve();
//    s.setBigPrimes();
//    s.countBigPrimes();
//    s.printSieveArray();
//    return 0;
//}

int main() {
    long x = 30000000;
    QuadrantSieve s(x);
    s.setMemberVariables();
    s.setSieveArray();
    s.printSieveArrayInfo();
    s.setSmallPrimes();
    s.sieve();
    s.setBigPrimes();
    s.countBigPrimes();
    return 0;
}