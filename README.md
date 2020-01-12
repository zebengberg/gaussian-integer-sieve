# gaussian-prime-sieve

>A program to generate primes in the Gaussian integers using the Sieve of Eratosthenes.


## Table of Contents

- [Gaussian Integers](#gaussian-integers)
- [Install](#install)
- [Command line usage](#command-line-usage)
- [Python bindings](#python-bindings)
- [Algorithm](#algorithm)
- [Implementation](#implementation)
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
 
Contact me if you would like better support for compiling this project on linux. I could readily use the `libstdc++` instead of the default macOS standard library.

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
![Segmented Block](/assets/block.png)
```Python
# Plot the Gaussian primes in the first quadrant with norm up to 10000. Use the save flag to write the image to disk.
>>> p = gprimes(10000)
>>> p.visualize(save=True)
```
![First Quadrant](/assets/first_quadrant.png)
```python
# Plot the Gaussian primes as well as their associates with norm up to 1000.
>>> p = gprimes(1000)
>>> p.visualize(full_disk=True)
```
![Full Plane](/assets/full_plane.png)


## Algorithm

With prime generating sieves, we often work with a *sieve array* A and a set of *small primes* P. For each prime p in P,the sieve proceeds by crossing off multiples of p in the sieve array A. To implement this on a computer, we define an array of booleans indexed by elements of A. The state of each boolean tracks if the element has yet been crossed off by some prime p in P.

In broad strokes, sieving in Z[i] involves the following steps.

1. Initialize a 2-dimensional array A with index (a, b) corresponding to the Gaussian integer a + bi.  In this way, we can view this array as representing a finite subset of the Gaussian integers.  Multiplication by i gives a 4-fold rotational symmetry in the complex plane; complex conjugation (in conjunction with multiplication by i) gives a reflectional symmetry about the line y = x in the complex plane. Such symmetries allow us to reduce the size of the sieving array if it contains redundancies encapsulated by these symmetries.

2. Generate a set of *small primes* P either through another sieving alogorithm or by using the sieve array to find prime elements of small norm. The set P should contain all primes with norm up to sqrt(B), where B is the norm of the largest Gaussian integer in A. In either case, for each prime a + bi in P, leave the sieve array A along at a + bi but mark off every Gaussian integer multiple of a + bi. Continue this process for each prime in P.

3. The unmarked entries that remain in the array A are exactly the Gaussian primes in A.


[See here](assets/algorithm.pdf) for a discussion of different approaches to generating primes in Z[i] and some motivation for our choice of algorithm.


### Crossing off multiples

Given a Gaussian integer a + bi and a sieve array A (A represents a sector or rectangle), we need to cross off multiples of a + bi from A. Said differently, we must cross off any sieve array element with index correspoinding to (a + bi)(c + di) for c, d in Z. Expanding this, we seek indices of the form (u, v) where
                u = ac - bd,  and  v = ad + bc.
This simply requires a double for-loop (once over c and once over d). Moreover, we can avoid costly multiplications by adding or subtracting d's or c's from the current index u and v when iterating over these for-loops. In the language of abstract algebra, we are crossing off elements of the Z-module generated by (a, b) and (-b, a).

### Donut Sieve

With the rational integers Z, the [wheel sieve](https://en.wikipedia.org/wiki/Wheel_factorization) allows for a more efficient way to manage and step through the sieve array. Primes other than 2, 3, and 5 must lie in one of the 8 distinct residue classes 1, 7, 11, 13, 17, 19, 23, 29 mod 30. In light of this observation, the wheel sieve (with a wheel size of 30) reduces both the memory and complexity of sieving by only tracking these 8 distinguished residue classes. If implemented carefully, a wheel sieve will not add significant overhead to the actual sieving process, and it will reduce to the total size of the sieve array to 8 / 30 = 27% of the original size.

We extend this idea of a wheel sieve to the Gaussian integers. Instead of a wheel corresponding to Z/30Z rolling over the integers, we consider the finite set of remainders in Z[i]/10Z[i]. Mathematically, this quotient is a **donut** or torus. Just as the primes 2, 3, 5 in a wheel sieve gave rise to a mod-30 wheel, the mod-10 donut arises from the primes 1 + i, 1 + 2i and 2 + i. The first of these three primes sits above 2, and the second both sit above 5.
- Every Gaussian integer z satisfying z = 0 (mod 1 + i) should be removed from the sieve array. The prime 1 + i has norm 2 and so one out of every two Gaussian integers will be divisible by 1 + i. In this way, including 1 + i in the donut will reduce the sieve array storage and memory by 50%.
- We deal with the primes 1 + 2i and 2 + i together (by the CRT, this is equivalent to working in Z[i]/5Z[i]). Every Gaussian integer satisfying z = 0 (mod 1 + 2i) **or** z = 0 (mod 2 + i). By the inclusion-exclusion principle and the Chinese remainder theorem, there should be exactly 9 solutions to these equations in Z[i]/5Z[i]. The Gaussian integer 5 has norm 25, and so 9 out of every 25 Gaussian integers will be removed from the sieve array. This provides a savings of 9/25 = 36%.
- Combining both of these local conditions, our donut allows us to only keep 32 of the 100 residue classes in Z[i]/10Z[i].

Just as in the wheel sieve, adding more primes to the wheel will continue to make the sieve array more efficient at the cost of increasing the combinatorial complexity of the donut itself. Although we stopped with a donut of size 10 for other reasons (see [implementation](c++-implementation), additional primes to consider include:
- using 3 in the donut would give a savings of 1/9 = 11% in sieve array size,
- using the primes 3 + 2i and 3 + 2i, the primes dividing 13, would give a savings of 25/169 = 15% in sieve array size.

In `DonutSieve.cpp` we implement the basic structures of this donut sieve. In particular, just as one hard codes gaps between wheel spokes in wheel sieve, the `DonutSieve` class contains certain variable and methods for access and dealing with the donut. As a simple example from this class, the array below enumerates the 32 residue classes we keep in the mod-10 donut.
```
{{-, 0, -, 1, -, -, -, 2, -, 3},
 {4, -, -, -, 5, -, 6, -, -, -},
 {-, -, -, 7, -, 8, -, 9, -, -},
 {10, -, 11, -, -, -, -, -, 12, -},
 {-, 13, -, -, -, 14, -, -, -, 15},
 {-, -, 16, -, 17, -, 18, -, 19, -},
 {-, 20, -, -, -, 21, -, -, -, 22},
 {23, -, 24, -, -, -, -, -, 25, -},
 {-, -, -, 26, -, 27, -, 28, -, -},
 {29, -, -, -, 30, -, 31, -, -, -}}
```

### Segmented Sieve

In a segmented sieve, we break up the entire sieving array into smaller subsets called **segments**. With these segments in place, we perform sieving (crossing off multiples of small primes) segment by segment. There are advantages and disadvantages to this segmented approach.
- The main disadvantage will be some additional overhead. If a python template for segmented sieving is:
```python
for segment in sieve_array:
    for p in small_primes:
        cross_off_multiples(segment, p)
```
then the inner for-loop will be smaller if the entire sieve array was used but it will also be repeated many times.
- The main advantage is that each segment can fit more readily into computer hardware such as the L2 cache. This will allow the inner for-loop in the snippet above to run more quickly because the sieve array segmented is cheap to access.

In this project, segmentation can be done by calling methods in the `BlockSieve` class.




## C++ Implementation

The aforementioned algorithm is implemented in a C++ library. `BaseSieve.cpp` defines an abstract base class with some of the basic sieving methods. Examples of classes derived from this include `QuadrantSieve`, `OctantSieve`, `DonutSieve`, and `SegmentedSieve`. Each derived class has its own method for initiating and accessing the sieve array. See the [usage examples](#command-line-usage) for various text representations of these sieve arrays.


Results of this repository can be verified with [http://oeis.org/A091100](http://oeis.org/A091100) to check for accuracy.

In the `OctantSieve` class, the 2-dimensional sieve array A is stored as a `vector<vector<bool>>` container. This C++ object is memory-efficient: each boolean is stored as a single bit in memory (this in itself gives a 224-fold savings over python 3.7 which requires 28 bytes of memory to store a boolean value).

In implementing the donut sieve, each 10 x 10 block of Gaussian integers corresponds to a full *donut roll*. The donut sieve requires us to track 32 residue classes for every 10 x 10 block of Gaussian integers. Said differently, every 10 x 10 block of Gaussian integers requires 32 bits of information to store their current state in the sieve process. Conveniently, a C++ `int` typically also requires 32 bits (4 bytes) of memory space. In this way, in donut-based classes such as `DonutSieve`, our sieve array is a `vector<vector<int>>` container in C++.





## Applications

This library can be used to generated data related to two unsolved problems in number theory.
- Prime number races in sectors
- Gaussian moat



## License

This project is released under the [MIT license](https://opensource.org/licenses/MIT).