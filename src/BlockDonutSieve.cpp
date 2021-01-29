/* Perform sieving in the box defined by [x, x + dx) x [y, y + dy). Must have
 * x > 0, y >= 0, dx > 0, dy > 0. Using a donut sieve approach; all of x, y, dx,
 * and dy must be multiples of 10.
 */

#include <iostream>
#include <cmath>
#include "OctantDonutSieve.hpp"
#include "BlockDonutSieve.hpp"
using namespace std;

// Using an initializer list
BlockDonutSieve::BlockDonutSieve(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy, bool verbose)
    : SieveTemplate<uint32_t>(pow((uint64_t)(x + dx - 1), 2) + pow((uint64_t)(y + dy - 1), 2), verbose), x(x) // x-coordinate of lower left-hand corner of segment block
      ,
      y(y) // y-coordinate of lower left-hand corner of segment block
      ,
      dx(dx) // horizontal side length of segment block
      ,
      dy(dy) // vertical side length of block

      // For a fixed c, d will iterate by d += 2 or d += 4, similar to the classic
      // wheel sieve. Got these iteration patterns in printDonut(), then pasted them
      // here so they can be hard coded in constructor.

      // donut array tracks how to move from one d to another during crossOffMultiples()
      // note: rows and columns reversed from usual complex plane orientation -- don't get confused.
      ,
      gapDonut{{0, 2, 0, 4, 0, 0, 0, 2, 0, 2},
               {4, 0, 0, 0, 2, 0, 4, 0, 0, 0},
               {0, 0, 0, 2, 0, 2, 0, 6, 0, 0},
               {2, 0, 6, 0, 0, 0, 0, 0, 2, 0},
               {0, 4, 0, 0, 0, 4, 0, 0, 0, 2},
               {0, 0, 2, 0, 2, 0, 2, 0, 4, 0},
               {0, 4, 0, 0, 0, 4, 0, 0, 0, 2},
               {2, 0, 6, 0, 0, 0, 0, 0, 2, 0},
               {0, 0, 0, 2, 0, 2, 0, 6, 0, 0},
               {4, 0, 0, 0, 2, 0, 4, 0, 0, 0}}

      ,
      bitDonut{{99, 0, 99, 1, 99, 99, 99, 2, 99, 3},
               {4, 99, 99, 99, 5, 99, 6, 99, 99, 99},
               {99, 99, 99, 7, 99, 8, 99, 9, 99, 99},
               {10, 99, 11, 99, 99, 99, 99, 99, 12, 99},
               {99, 13, 99, 99, 99, 14, 99, 99, 99, 15},
               {99, 99, 16, 99, 17, 99, 18, 99, 19, 99},
               {99, 20, 99, 99, 99, 21, 99, 99, 99, 22},
               {23, 99, 24, 99, 99, 99, 99, 99, 25, 99},
               {99, 99, 99, 26, 99, 27, 99, 28, 99, 99},
               {29, 99, 99, 99, 30, 99, 31, 99, 99, 99}}

      ,
      realPartDecompress{0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9}, imagPartDecompress{1, 3, 7, 9, 0, 4, 6, 3, 5, 7, 0, 2, 8, 1, 5, 9, 2, 4, 6, 8, 1, 5, 9, 0, 2, 8, 3, 5, 7, 0, 4, 6}
{
  if ((x % 10) || (y % 10) || (dx % 10) || (dy % 10))
  {
    cerr << "All in parameters should be multiples of 10 to use this donut method." << endl;
    throw;
  }
}

// Calling the fastest sieve method available to generate small primes.
void BlockDonutSieve::setSmallPrimes()
{
  if (verbose)
  {
    cerr << "Calling the OctantDonutSieve to generate smallPrimes..." << endl;
  }
  OctantDonutSieve s((uint64_t)isqrt(maxNorm), false);
  s.run();
  smallPrimes = s.getBigPrimes();
}

// Sieve array holds indices corresponding to gints with x <= a < x + dx and
// y <= b < y + dy. Since we are using a donut approach, the size of the sieve
// array is dx / 10 by dy / 10.
void BlockDonutSieve::setSieveArray()
{
  if (verbose)
  {
    cerr << "Building donut sieve array..." << endl;
  }
  // trick to take integer division ceiling instead of floor
  // each donut block must start at a multiple of 10
  for (uint64_t i = 0; i < dx / 10; i++)
  {
    vector<uint32_t> column(dy / 10, (uint32_t)pow(2, 32) - 1);
    sieveArray.push_back(column);
  }
  if ((x == 0) && (y == 0))
  {
    setFalse(1, 0); // Crossing off 1
    setFalse(0, 1); // Crossing off i
  }
  if (verbose)
  {
    printSieveArrayInfo();
  }
}

// Cross off multiples of the gint g = a + bi within the sieveArray.
// Let a + bi be the gint and c + di be the co-factor of the multiple we seek.
// We must over integers c and d such that the product (a + bi)(c + di) is
// contained in the segment block. This gives inequalities x <= ac - bd < x + dx
// and y <= ad + bc < y + dy. Solve these for c to get
// (ax + by) / (a^2 + b^2) <= c <= (a(x + dx - 1) + b(y + dy - 1)) / (a^2 + b^2).
// Then d is defined by some max and min conditions.
void BlockDonutSieve::crossOffMultiples(gint g)
{
  if (g.norm() <= 5)
  {
    return;
  } // exit early if g is above 2 or 5.
  // Convert everything to unsigned long type.
  int64_t a = g.a;
  int64_t b = g.b;
  int64_t N = g.norm();
  int64_t c, cUpper;
  // Trick: to get ceiling of integer division n / m, call (n + m - 1) / m.
  // This only works when both m and n > 0.
  if (g.b)
  {
    c = (a * x + b * y + N - 1) / N; // using trick
    cUpper = (a * (x + dx - 1) + b * (y + dy - 1)) / N;
  }
  else
  {
    c = (x + a - 1) / a; // using trick
    cUpper = (x + dx - 1) / a;
  }
  for (; c <= cUpper; c++)
  {
    int64_t d, dUpper; // d isn't necessarily positive
    if (b)
    {
      // trick doesn't work because of sign; instead using ceil.
      double d1 = double(a * c - x - dx + 1) / double(b);
      double d2 = double(y - b * c) / double(a);
      d = int64_t(ceil(max(d1, d2)));
      // Reusing d1 and d2
      d1 = double(a * c - x) / double(b);
      d2 = double(y + dy - 1 - b * c) / double(a);
      dUpper = int64_t(floor(min(d1, d2)));
    }
    else
    {
      d = (y + a - 1) / a; // using trick
      dUpper = (y + dy - 1) / a;
    }
    // Annoying indexing because d can be negative; need to use mod function
    // defined in the BaseSieve.cpp file.
    int32_t jump = gapDonut[c % 10][mod(d, 10)];
    // Now trying to figure out where to start d so that c + di is coprime to 10.
    while (jump == 0)
    {
      d++;
      jump = gapDonut[c % 10][mod(d, 10)];
    }
    int64_t u = a * c - b * d - x; // u = ac - bd - x
    int64_t v = b * c + a * d - y; // v = bc + ad - y
    while (d <= dUpper)
    {
      setFalse(u, v);
      jump = gapDonut[c % 10][mod(d, 10)];
      d += jump;
      u -= jump * b;
      v += jump * a;
    }
  }
  if ((x <= a) && (a < x + dx) && (y <= b) && (b < y + dy))
  {
    // crossed this off; need to re-mark it as prime
    setTrue(a - x, b - y);
  }
  if ((x <= b) && (b < x + dx) && (y <= a) && (a < y + dy))
  {
    // crossed this off; need to re-mark it as prime
    setTrue(b - x, a - y);
  }
}

// Set the bit in the sieveArray to false corresponding to the gint u + vi
void BlockDonutSieve::setFalse(uint32_t u, uint32_t v)
{
  uint32_t bit = bitDonut[u % 10][v % 10];
  // clearing the bit; 1u is unsigned int
  sieveArray[u / 10][v / 10] &= ~(1u << bit);
}

// Set the bit in the sieveArray to true corresponding to the gint u + vi
void BlockDonutSieve::setTrue(uint32_t u, uint32_t v)
{
  unsigned char bit = bitDonut[u % 10][v % 10];
  // setting the bit to 1; 1u is unsigned int
  sieveArray[u / 10][v / 10] |= (1u << bit);
}

// Combing through the sieve array to find primes after sieving.
void BlockDonutSieve::setBigPrimes()
{
  if (verbose)
  {
    cerr << "Gathering primes after sieve..." << endl;
  }
  // Putting in primes dividing 10.
  if ((x < 10) && (y < 10))
  {
    bigPrimes.emplace_back(1, 1);
    bigPrimes.emplace_back(2, 1);
    bigPrimes.emplace_back(1, 2);
  }

  for (uint32_t a = 0; a < dx / 10; a++)
  {
    for (uint32_t b = 0; b < dy / 10; b++)
    {
      for (unsigned char bit = 0; bit < 32; bit++)
      {
        if ((sieveArray[a][b] >> bit) & 1u)
        {
          gint g(x + 10 * a + realPartDecompress[bit],
                 y + 10 * b + imagPartDecompress[bit]);
          bigPrimes.push_back(g);
        }
      }
    }
  }
  if (verbose)
  {
    cerr << "Done gathering." << endl;
  }
}

uint64_t BlockDonutSieve::getCountBigPrimes()
{
  if (verbose)
  {
    cerr << "Counting primes after sieve..." << endl;
  }
  uint64_t count = 0;
  // Putting in primes dividing 10.
  if ((x == 0) && (y == 0))
  {
    count = 3;
  }
  for (uint32_t a = 0; a < dx / 10; a++)
  {
    for (uint32_t b = 0; b < dy / 10; b++)
    {
      for (unsigned char bit = 0; bit < 32; bit++)
      {
        if ((sieveArray[a][b] >> bit) & 1u)
        {
          count++;
        }
      }
    }
  }
  if (verbose)
  {
    cerr << "Total number of primes: " << count << "\n"
         << endl;
  }
  return count;
}