#include <iostream>
#include <cmath>
#include "gint.hpp"

using namespace std;

class OctantSieve{
private:
    long x;
    vector<vector<bool>> sieveArray;
    vector<Gint> primes;  // to hold primes after sieving

public:
    explicit OctantSieve(long x) { this-> x = x; };  // constructor
    void initializeSieveArray();
    void crossOffMultiplesSkew(Gint);
    void crossOffMultiplesRect(Gint);
    void getPrimes();
    void printPrimes();
    void writePrimesToFile();
    void printProgress(Gint);  // update after each sieved prime
};

void OctantSieve::initializeSieveArray() {
    // sieveArray holds values for Gint's with a, b >= 0, a >= b, and a^2 + b^2 <= x.
    cout << "Building sieve array..." << endl;
    for (long a = 0; a <= isqrt(x); a++) {
        // Calculating the intersection of circle a^2 + b^2 <= x and the line a = b.
        long intersection = long(sqrt(double(x) / 2.0));
        long b = a <= intersection ? a + 1 : isqrt(x - a * a) + 1;
        vector<bool> column((unsigned long)b, true);  // Create a vector of size b with all values true.
        sieveArray.push_back(column);
    }
    sieveArray[0][0] = false;  // 0 is not prime
    sieveArray[1][0] = false;  // 1 is not prime
}

int main() {
    OctantSieve os(1000);
    os.initializeSieveArray();
}