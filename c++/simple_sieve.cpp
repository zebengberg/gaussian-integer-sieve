#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct point {
    int a, b, norm;
    point(int a, int b) {this->a = a; this->b = b; this->norm = a * a + b * b;}
};

// Read small primes to be used for sieving from existing file.
// Getting all primes from file small_primes.txt with norm up to x.
vector<point> readPrimesFromFile(int x) {
    // If x is too small, then we are not actually taking any primes.
    if (x < 2) {
        cerr << "No primes read and no sieving needed.";
        exit(1);
    }
    // small_primes.txt contains a list of Gaussian primes sorted by norm.
    ifstream f;
    f.open("/Users/robotics/drive/git/gaussian-prime-sieve/data/small_primes.txt");
    if (!f) {  //TODO: try to make this a relative path
        cerr << "Unable to open file small_primes.txt";
        exit(1);   // call system to stop
    } else {
        cout << "Reading small primes from file." << endl;
    }

    int a, b;
    f >> a >> b;  // token-based parsing
    vector<point> smallPrimes;  // the list of small primes to be populated
    while (a * a + b * b <= x) {
        point p(a, b);  // creating prime struct
        smallPrimes.push_back(p);  // putting it into vector P

        f >> a >> b;  // reading next pair a, b from file f
        // If we've gotten to the end of the file, then we don't have enough precomputed primes.
        if (f.eof()) {
            cerr << "Not enough primes in small_primes.txt";
            exit(1);
        }
    }
    f.close();
    return smallPrimes;
}

// Integer square root.
int isqrt(int n) {
    int x, y;
    x = n;
    y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

class SieveArray{
        int x;
        vector<vector<bool>> sieveArray;
        vector<point> primes;
    public:
        void initializeSieveArray(int);
        void crossOffMultiples(point);
        void getPrimes();
        void printPrimes();
};

void SieveArray::initializeSieveArray(int x) {
    // Lots of controversy about vector<bool>.
    // See https://stackoverflow.com/questions/17794569/why-is-vectorbool-not-a-stl-container
    // Each boolean value is stored as a single bit but pointers are not available.
    this->x = x;
    for (int a = 0; a <= isqrt(x); a++) {
        int b = isqrt(x - a * a) + 1;
        vector<bool> row(b, true);  // Create a vector of size b with all values true.
        sieveArray.push_back(row);
    }
}

void SieveArray::crossOffMultiples(point p) {
    int outer_u = 0;
    int outer_v = 0;
    for (int c = 0; c <= isqrt(x); c++) {
        int u = outer_u;
        int v = outer_v;
        for (int d = 0; d <= isqrt(x / p.norm - c * c); d++) {
            if (u > 0) {
                sieveArray[u][v] = false;
            } else {
                sieveArray[v][-u] = false;
            }
            u -= p.b;
            v += p.a;
        }
        outer_u += p.a;
        outer_v += p.b;
    }
}


void SieveArray::getPrimes() {
    for (int a = 0; a <= isqrt(x); a++) {
        for (int b = 0; b <= isqrt(x - a * a); b++) {
            if (sieveArray[a][b]) {
                point p(a, b);
                primes.push_back(p);
            }
        }
    }
}

void SieveArray::printPrimes() {
    for (point p : primes) {
        cout << p.a << "  " << p.b << "  " << p.norm << endl;
    }
    cout << "Total number of primes: " << primes.size() << endl;
}


int main(){
    SieveArray sieveArray;
    sieveArray.initializeSieveArray(10000);
    vector<point> smallPrimes = readPrimesFromFile(100);
    for (point p : smallPrimes) {
        sieveArray.crossOffMultiples(p);
    }
    sieveArray.getPrimes();
    sieveArray.printPrimes();
}

// TODO: something wrong with small primes read from file. Should get 1232 for count with norm <= 10000,