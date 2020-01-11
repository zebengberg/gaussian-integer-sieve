# gaussian-prime-sieve

>A program to generate primes in the Gaussian integers using the Sieve of Eratosthenes.


## Table of Contents

- [Gaussian Integers](#gints)
- [Install](#install)
- [Command line usage](#usage)
- [Python bindings](#bindings)
- [Sieving](#sieving)
- [Algorithm](#algorithm)
- [Applications](#appplications)
- [License](#license)


## Gaussian Integers

The **Gaussian integers** are complex numbers of the form a + bi for integers a and b. Here i is the **imaginary unit** with the property that i * i = -1. The set of Gaussian integers, denoted by Z[i], forms a mathematical structure that enjoys many of the same properties as those of the rational integers, denoted by Z. Gaussian integers can be factored and the notion of primality is well-defined. The **norm** of a Gaussian integer gives a measure of its size; the analogy in the rational integers is absolute value. With these tools in hand, sieving ideas can be used to compute primes in this extension of our usual integer system.

This project is an implementation of the *Sieve of Eratosthenes* to generate primes in the Gaussian integers. While the Gaussian integers present some difficulties that mostly nonexistent in the context of the rational integers, nearly all of the basic computational sieve ideas can be extended to Z[i]. In particular, both the wheel sieve and the segmented sieve can be extended to the Gaussian integers.

The goal of this project is to build an efficient implementation of a prime generating sieve in the Gaussian integers in both Python and C++, and to provide a library for deeper exploration. The [*primesieve*](https://github.com/kimwalisch/primesieve) project implements the current state of the art prime generator for rational integers. *primesieve* includes both algorithmic optimizations as well as hardware considerations needed in order to improve runtime performance. In this repository, I hope to extend some of the more accessible algorithmic and hardware-optimization ideas from *primesieve*.


## Install

This entire repository can be cloned or downloaded for use. To build the executable on a macOS, cd into the
 `gaussian-integer-sieve` directory and run
```shell script
$ make
```
This requires a C++ compiler supporting both C++11 and the `libc++` library as well as `make`. On macOS, they can be
 installed with `xcode-select --install`.

 The Cython bindings can be built from source with
 ```shell script
$ python setup.py build_ext --inplace
```
from within the `gaussian-integer-sieve/python` directory. Building this module requires Python 3.x, Cython, and the
 aforementioned C++ compiler. After the module is compiled, it can be imported into a python console or used in a
  python script.


## Command line usage
Command line options include:
```
Usage: ./gintsieve x [y z] [option1] [option2] ...
Generate Gaussian primes with norm up to x using sieving methods.
    x                   The norm-bound of the generated primes
    y                   Coordinates (x, y) of SW-corner of block in segmented sieve.
    z                   Side length block in segmented sieve.

Options:
    -h, --help          Print this help message.
    -v, --verbose       Display sieving progress.
    -p, --printprimes   Print the real and imag part of primes found by the sieve.
    -w, --write         Write primes to csv file in working directory.
    -a, --printarray    Print a text representation of the sieve array.
    -c, --count         Count the number of generated primes and exit program.
    -q, --quadrant      Sieve array consists of Gaussian integers in the first quadrant.
    -o, --octant        Sieve array consists of Gaussian integers in the first octant.
    -d, --donut         Sieve array consists of Gaussian integers in first octant
                        coprime to 2 and 5.
    -s, --segmented     Sieve array consists of Gaussian integers of form a + bi with
                        x <= a < x + z and y <= b < y + z.
```
For example, to get the real and imaginary parts of the Gaussian primes with norm up to 60 sorted by norm, run:
```shell script
$ ./gintsieve 60
1 1
2 1
1 2
3 0
3 2
2 3
4 1
1 4
5 2
2 5
6 1
1 6
5 4
4 5
7 0
7 2
2 7
Total number of primes, including associates: 68
```
Many of the sieving algorithms can be visualized by using the `-a` option to print the sieve array after the sieving
 is finished. For example, we can see a hexidecimal encoding of the donut sieving with:
 ```shell script
$ ./gintsieve 15000 -d -a
                                                                        bf6edbef ffffffec                            
                                                               434cdbff e01dbe2c db48e06c ffffffff                   
                                                      0aedf3ff c0120c78 88730b2e 0aae83e5 fffff145                   
                                             6fdcd7ff 1f08ddf0 eed55590 d18b13c7 35e62db2 c3496c47 fffffffc          
                                    f65cffef 70113512 369847a3 36968225 e82622d5 1149408b 9490ea99 fffede51          
                           43eedbff f8cba42e 51413f89 082e18b5 88face6d 2a0c2722 5316d330 fb516844 ff0946a4          
                  0bffdfef 0e47b6da 25f0b9da f9245714 87cac0ea 8894d858 8f4c51d5 9602fc1a 6914b816 792d101d fffffffe 
         f3ceffef ac979b9b be4c61d3 94f6cea5 64332c6d 660ea2b6 1196506d 40fdb31f 2a1e034f 35c0ac56 08b210f2 ffffffd0 
5fffffee bb73ba7f c9767527 1739d9f0 c3a99523 bf5e23a4 199f8080 f06afb52 554124ac 4f237261 48e00f02 c5530187 fffffc25 
```
To see the primality of a 15 x 15 block of Gaussian integers way out the in complex plane, we run:
```shell script
$ ./gintsieve 100000000 300000000 15 -a -s
---------------
---------------
-------------*-
--------------*
---*-------*---
---------------
---------------
---------------
---------*-----
------*--------
---------------
---------------
---*---------*-
---------------
-*-------------
```

## Python bindings

Once the python module is built, we can import it into python and use it as in the following examples.


```Python
>>>  from gintsieve import *

# Generate a list of tuples representing Gaussian primes with norm up to 50 sorted by norm.
>>> gprimes(40)
[(1, 1), (2, 1), (1, 2), (3, 0), (3, 2), (2, 3), (4, 1), (1, 4), (5, 2), (2, 5), (6, 1), (1, 6)]

# Instead of returning a list, we can count the Gaussian primes.
>>> gprimes_count(3141592653)
# Calling the QuadrantSieve to generate smallPrimes...
# Building sieve array...
# Sieve array approximate memory use: 49MB.
# Starting to sieve...
# Done sieving. Total time for sieving: 15.499 seconds.
# Counting primes after sieve...
# Done with count.
# Total number of primes, including associates: 603726004
603726004

# The function gprimes actually returns an object that can be plotted. Below we get a
# block of Gaussian primes generated with a segmented sieve.
>>> p = gprimes(123456, 67890, 100)
>>> p.visualize()
```
![Segmented Block](/images/block.png)
```Python
# Plot the Gaussian primes in the first quadrant with norm up to 10000. Use the save flag to write the image to disk.
>>> p = gprimes(10000)
>>> p.visualize(save=True)
```
![First Quadrant](/images/first_quadrant.png)
```python
# Plot the Gaussian primes as well as their associates with norm up to 1000.
>>> p = gprimes(1000)
>>> p.visualize(full_disk=True)
```
![Full Plane](/images/full_plane.png)


## Sieving




With prime generating sieves, we often work with a *sieve array* A and a set of primes P. For each prime p in P, the
 sieve proceeds by crossing off multiples of p in the sieve array A. To implement this on a computer, we define an
  array of booleans indexed by elements of A. The state of each boolean tracks if the element has yet been crossed
   off by some prime p in P.

Working within the Gaussian integers, our sieve array A is a 2-dimensional object indexed by Gaussian integers. In
 native python, this can be implemented with a list of lists. We can also use a 2-dimensional *numpy* array. In C++, we can use the built-in array, or vector, or another class of container.

In python 3.7, the boolean `True` takes up 28 bytes of memory storage. The vast majority of this is overhead; in theory, a boolean value should only occupy a single bit of space. Even in C++, a single boolean value requires a full byte of memory. One issue involved with storing a boolean value into a smaller storage space is that modern computers do not allow pointers to individual bits in memory. In C++, there are various ways around this issue such as the `bitset` object.

A different way to resolve this issue is to include some of the ideas of the wheel sieve. Primes other than 2, 3, and 5 must lie in one of the 8 distinct residue classes 1, 7, 11, 13, 17, 19, 23, 29 mod 30. In light of this observation, the wheel sieve reduces both the memory and complexity of sieving by only tracking these 8 distinguished residue classes. If implemented carefully, a wheel sieve should not add significant overhead to the actual sieving process. Moreover, this approach dovetails nicely with modern computer architecture: a single byte can track the state of 30 consecutive integers from an integer sieve array, allowing for a 30-fold reduction in memory usage (in a C++ implementation). Extending this wheel sieve to the Gaussian integers, both the memory and complexity of the sieve can be greatly reduced.

In `simple_sieve.py` and `simple_sieve.cpp`, I use native booleans to track the state of each Gaussian integer in the sieve array. This will create significant overhead in both python (up to 28 * 8 = 224-fold waste) and C++ (8-fold).  In our donut sieve (TODO: implement this), I use the aforementioned ideas to reduce and align the storage requirements of the sieve to modern hardware architecture.

Results of this repository can be verified with [http://oeis.org/A091100](http://oeis.org/A091100) to check for accuracy.


## Algorithms

There are several immediate methods that can be used to find all prime numbers in **Z[i]** up to norm **x**. We describe
 three below.
 1. Using the sieve of Eratosthenes in **Z**, generate all prime numbers up to **x** as some array **A**. According to
  elementary number theory, a Gaussian integer **a + bi** with **a > 0** and **b >= 0** is prime if and only if **N(a
   + bi)** is a prime in **Z** or **b = 0** and **a** is a prime in **Z**. Using the sieve takes **O(x log log x
   )** steps, and checking each Gaussian integer for primality in this way takes **O(x)** steps. The entire runtime of
    this algorithm is **O(x log log x)**.
 2. 






The module `simple_sieve.py` can be used to generate Gaussian integer primes with norm up to x. This algorithm is summarized in the following steps.

1. Initialize a 2-dimensional array with index (a, b) corresponding to the Gaussian integer a + bi.  In this way, we can simply view this array as a finite subset of the Gaussian integers.  By restricting to the first quadrant of the complex plane, we take a and b to be non-negative integers with a^2 + b^2 <= x.  Mark off the array at the Gaussian integers 0, 1, and i.  These Gaussian integers are neither prime nor composite.
2. Locate the unmarked Gaussian integer a + bi; with smallest norm.  Necessarily a + bi is prime.  Leave the array empty at a + bi but mark off every Gaussian-integer-multiple of a + bi.
3. Loop over the process in step 2 by marking off multiples of other previously unmarked Gaussian primes.  Continue looping until every unmarked entry with norm less than sqrt x has been used for sieving.
4. The unmarked entries that remain in the array are exactly the Gaussian integer primes with norm up to x.

Compared with the usual sieve of Eratosthenes in the rational integers, step 2 is somewhat more involved in Z[i].  There are several naive approaches to this.

- Crossing off multiples of the Gaussian integer a + bi is equivalent to finding elements of the Z[i]-lattice generated by a + bi. In turn, this is equivalent to finding elements of the Z-lattice generated by the points (a, b) and (-b, a). This approach is implemented in the function `cross_off_multiples()`.

- In the case that a + bi is a degree 1 prime lying over p, multiples of a + bi form a subgroup of order p in the quotient group Z[i] / pZ[i]. Starting with representatives of this subgroup, we can then translate these representatives in the real and imaginary directions by multiples of p. Because this subgroup has order p, it is cyclic, and representative can be attained by taking Z-multiples of a + bi (instead of Z[i]-multiples). This approach is implemented in the function `cross_off_multiples2()`.

These two approaches have trade-offs. The first approach can be implemented with two nested for-loops, whereas three such for-loops are needed in the second approach. In the first approach, translating through the sieve array is not done parallel to the coordinate axes, and so more care must be exercised to determine when the index is out of bounds. In the second approach, translation is always occuring within a particular row or column of the sieve array. Generally, the second approach is faster for small primes whereas the first approach is faster for large primes. The choice is which of these two functions to use is not currently optimized.



The program `simple_sieve.cpp` takes a similar approach to the implementation in python. The only difference here is that `simple_sieve.cpp` reads a file containing a list of small primes used to sieve rather than generating such primes on the fly (as is done in `simple_sieve.py`). This C++ implementation requires a list of primes with norm up to sqrt(x) so that it can sieve to compute all primes with norm up to x.

## Applications



## License

This project is released under the [MIT license](https://opensource.org/licenses/MIT).