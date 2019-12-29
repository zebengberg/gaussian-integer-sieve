#ifndef C___GINT_HPP
#define C___GINT_HPP

#include <vector>
using namespace std;

// Gaussian integer class.
class Gint {
private:
    long real, imag;
public:
    Gint(long a, long b) {  // constructor
        real = a;
        imag = b;
    }
    long a() { return real; }
    long b() { return imag; }
    long norm() { return real * real + imag * imag; }
};

// Functions defined in QuadrantSieve.cpp.
vector<Gint> readPrimesFromFile(long);
long isqrt(long);


#endif //C___GINT_HPP
