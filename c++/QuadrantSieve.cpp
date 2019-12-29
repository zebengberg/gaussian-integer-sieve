#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>

using namespace std;

// Gaussian integer class.
class Gint {
private:
    int real, imag;
public:
    Gint(int a, int b) {  // constructor
        real = a;
        imag = b;
    }
    int a() { return real; }
    int b() { return imag; }
    int norm() { return real * real + imag * imag; }
};

// Read small primes to be used for sieving from existing file.
// Getting all primes from file small_primes.txt with norm up to x.
vector<Gint> readPrimesFromFile(int x) {
    // If x is too small, then we are not actually taking any primes.
    if (x < 2) {
        cerr << "No primes read and no sieving needed.";
        exit(1);
    }
    // small_primes.txt contains a list of Gaussian primes sorted by norm.
    ifstream f;
    f.open("../data/small_primes.txt");
    if (!f) {
        cerr << "Unable to open file small_primes.txt";
        exit(1);   // call system to stop
    } else {
        cout << "Reading small primes from file." << endl;
    }

    int a, b;
    f >> a >> b;  // token-based parsing
    vector<Gint> smallPrimes;  // the list of small primes to be populated
    while (a * a + b * b <= x) {
        Gint p(a, b);  // creating prime struct
        smallPrimes.push_back(p);  // putting it into vector P

        f >> a >> b;  // reading next pair a, b from file f
        // If we've gotten to the end of the file, then we don't have enough precomputed primes.
        if (f.eof()) {
            cerr << "Not enough primes in small_primes.txt";
            f.close();
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

class QuadrantSieveArray{
private:
    int x;
    double progress;
    double totalProgress;
    vector<vector<bool>> sieveArray;
    vector<Gint> primes;

public:
    QuadrantSieveArray(int);  // constructor
    void initializeSieveArray();
    void crossOffMultiplesSkew(Gint);
    void crossOffMultiplesRect(Gint);
    void getPrimes();
    void printPrimes();
    void writePrimesToFile();
    void printProgress(Gint);  // update after each sieved prime
};

QuadrantSieveArray::QuadrantSieveArray(int x) {
    this->x = x;
    progress = 0;
    totalProgress = log(log(x)) - log(2.0);
}

void QuadrantSieveArray::initializeSieveArray() {
    // Lots of controversy about vector<bool>.
    // Each boolean value is stored as a single bit but pointers are not available.
    cout << "Building sieve array..." << endl;
    for (int a = 0; a <= isqrt(x); a++) {
        int b = isqrt(x - a * a) + 1;
        vector<bool> row(b, true);  // Create a vector of size b with all values true.
        sieveArray.push_back(row);
    }
    sieveArray[0][0] = false;  // 0 is not prime
    sieveArray[1][0] = false;  // 1 is not prime
    sieveArray[0][1] = false;  // i is not prime

    auto size = 3.1415 * x / 4;  // area of the quarter circle
    size /= 8;  // 8 bits per byte
    size /= 1000000000;  // convert to GB
    cout << "Sieve array approximate memory use: " << size  << "GB" << endl;
}

void QuadrantSieveArray::crossOffMultiplesSkew(Gint g) {
    for (int c = 1; c <= isqrt(x); c++) {
        int u = c * g.a();
        int v = c * g.b();
        for (int d = 0; d <= isqrt(x / g.norm() - c * c); d++) {
            if (u > 0) {
                sieveArray[u][v] = false;
            } else {  // multiply u + vi by -i to bring it back to the first quadrant
                sieveArray[v][-u] = false;
            }
            u -= g.b();
            v += g.a();
        }
    }
    sieveArray[g.a()][g.b()] = true;  // crossed this off; need to remark it as prime
}

void QuadrantSieveArray::crossOffMultiplesRect(Gint g) {
    int p, subgroupSize;
    if (g.b()) {  // degree 1 primes
        p = g.norm();
        subgroupSize = p;
    } else {  // degree 2 primes
        p = g.a();
        subgroupSize = 1;  // trivial subgroup
    }
    for (int i = 0; i < subgroupSize; i++) {
        // s + it is an element of the additive subgroup < a + bi > in Z[i] / pZ[i]
        int s = i * g.a() % p;
        int t = i * g.b() % p;
        for (int u = s; u <= isqrt(x); u += p) {
            for (int v = t; v <= isqrt(x - u * u); v += p) {
                sieveArray[u][v] = false;
            }
        }
    }
    sieveArray[g.a()][g.b()] = true;  // crossed this off; need to remark it as prime
}

void QuadrantSieveArray::getPrimes() {
    cout << "Gathering primes after sieve..." << endl;
    for (int a = 1; a <= isqrt(x); a++) {  // Want to stay off imaginary line
        for (int b = 0; b <= isqrt(x - a * a); b++) {
            if (sieveArray[a][b]) {
                Gint p(a, b);
                primes.push_back(p);
            }
        }
    }
}

void QuadrantSieveArray::printPrimes() {
    for (Gint p : primes) {
        cout << p.a() << "  " << p.b() << "  " << p.norm() << endl;
    }
    cout << "Total number of primes: " << primes.size() << endl;
}

void QuadrantSieveArray::writePrimesToFile() {
    ofstream f;
    f.open("../data/cpp_primes.csv");
    for (Gint p : primes) {
        f << p.a() << "," << p.b() << "," << p.norm() << endl;
    }
    f.close();
}

void QuadrantSieveArray::printProgress(Gint g) {
    progress += 1.0 / g.norm();
    int barPos = 70 * progress / totalProgress;
    cout << "[";
    for (int i = 0; i < barPos; i++) { cout << "|"; }
    for (int i = barPos; i <= 70; i++) { cout << " "; }
    cout << "]\r";
    cout.flush();
}


int main(int argc, const char* argv[]){
    if (argc != 2) {
        cerr << "Specify the sieve array norm bound." << endl;
        return 1;
    }
    int x = atoi(argv[1]);  // convert command line input to int
    QuadrantSieveArray sieveArray(x);
    sieveArray.initializeSieveArray();
    vector<Gint> smallPrimes = readPrimesFromFile(isqrt(x));
    cout << "Starting to sieve..." << endl;
    auto startTime = chrono::high_resolution_clock::now();
    for (Gint p : smallPrimes) {
        sieveArray.crossOffMultiplesRect(p);
        sieveArray.printProgress(p);
    }
    auto endTime = chrono::high_resolution_clock::now();
    auto totalTime = chrono::duration_cast<chrono::seconds>(endTime - startTime);
    cout << "Total time for sieving: " << totalTime.count() << " seconds" << endl;
    sieveArray.getPrimes();
    //sieveArray.printPrimes();
    //sieveArray.writePrimesToFile();
    return 0;
}
