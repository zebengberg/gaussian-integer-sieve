#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct point {
    int a, b;
    point(int a, int b) {this->a = a; this->b = b;}
};

// Read small primes to be used for sieving from existing file.
// Getting all primes from file small_primes.txt with norm up to x.
vector<point> read_primes_from_file(int x) {

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
    vector<point> P;  // the list of small primes to be populated
    while (a * a + b * b <= x) {
        point p(a, b);  // creating prime struct
        P.push_back(p);  // putting it into vector P

        f >> a >> b;  // reading next pair a, b from file f
        // If we've gotten to the end of the file, then we don't have enough precomputed primes.
        if (f.eof()) {
            cerr << "Not enough primes in small_primes.txt";
            exit(1);
        }
    }
    f.close();
    return P;
}


int main(){
    vector<point> P = read_primes_from_file(100);
}
