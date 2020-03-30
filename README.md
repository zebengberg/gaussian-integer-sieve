# gaussian-prime-sieve

>A program to generate primes in the Gaussian integers using the Sieve of Eratosthenes.


## Table of Contents

- [Gaussian Integers](#gaussian-integers)
- [Install](#install)
- [Command line usage](#command-line-usage)
- [Python bindings](#python-bindings)
- [Algorithm](#algorithm)
- [C++ Implementation](#c++-implementation)
- [Applications](#applications)
- [License](#license)


## Gaussian Integers

The **Gaussian integers** are complex numbers of the form a + bi for integers a and b. Here i is the **imaginary unit**, so i * i = -1. The set of Gaussian integers, denoted by Z[i], forms a mathematical structure that enjoys many of the same properties as those of the **rational integers** (the usual integers), denoted by Z. Gaussian integers can be factored and the notion of primality is well-defined. The **norm** of a Gaussian integer gives a measure of its size; the analogy in the rational integers is absolute value. With these tools in hand, sieving ideas can be used to compute primes in this extension of our usual integer system.

This project is an implementation of the *Sieve of Eratosthenes* to generate primes in the Gaussian integers. While the Gaussian integers present some difficulties that mostly nonexistent in the context of the rational integers, nearly all of the basic computational sieve ideas can be extended to Z[i]. In particular, both the wheel sieve and the segmented sieve can be extended to the Gaussian integers.

The goal of this project is to build an efficient C++ implementation of a prime generating sieve in the Gaussian integers, and to provide a library for deeper exploration. The [*primesieve*](https://github.com/kimwalisch/primesieve) project implements the current state of the art prime generator for rational integers. *primesieve* includes both algorithmic optimizations as well as hardware considerations needed in order to improve runtime performance. In this repository, I draw inspiration from *primesieve* and hope to extend some of its algorithms and hardware-optimizations to the realm of the Gaussian integers.


## Install

This entire repository can be cloned or downloaded for use. To build the executable on macOS, cd into the
 `gaussian-integer-sieve` directory and run
```shell script
$ make
```
This requires a C++ compiler supporting both C++11 and the `libc++` library as well as `make`. On macOS, these can be installed by running `xcode-select --install`.
 
Please contact me if you would like better support for compiling this project on linux. Linux systems generally use `libstdc++` instead of the default macOS library `libc++`.

 The [Python bindings](#python-bindings) can be built from source with
 ```shell script
$ python setup.py build_ext --inplace
```
from within the `gaussian-integer-sieve/python` directory. Building this module requires Python 3.x, Cython, and the aforementioned C++ compiler. After the module is compiled, it can be imported into a python console or used in a python script.


## Command line usage

All sieving algorithms can be called from the `gintsieve` executable.
Command line options include:
```
Usage: ./gintsieve x [y dx dy alpha beta] [option1] [option2] ...
Generate Gaussian primes with norm up to x using sieving methods.
    x                   The norm-bound of the generated primes
    y                   Coordinates (x, y) of SW-corner of array in block sieve.
    dx                  Horizontal side length in block sieve.
    dy                  Vertical side length in block sieve.
    alpha               Start angle in sector sieve.
    beta                Final angle in sector sieve.

Options:
    -h, --help          Print this help message.
    -v, --verbose       Display sieving progress.
    -p, --printprimes   Print the real and imag part of primes found by the sieve.
    -w, --write         Write primes to csv file in current directory.
    -a, --printarray    Print a text representation of the sieve array.
    -c, --count         Count the number of generated primes and exit program.

Optional sieve types:
    -o, --octant        Sieve array indexed by Gaussian integers in the first octant.
                        This is the default sieve method called.
    -s, --sector        Sieve array indexed by Gaussian integers in the sector
                        with start angle alpha and final angle beta.
    -b, --block         Sieve array indexed by Gaussian integers in the rectangle
                        defined by x <= real < x + dx and y <= imag < y + dy.
    -d, --donut         If a donut version of the sieve array exists, use it. In the
                        donut sieve, the sieve array consists of Gaussian integers
                        coprime to 2 and 5. This optional can be used with --octant
                        and --block, and is often significantly faster.
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
Total number of primes printed: 17
```
Many of the sieving algorithms can be visualized by using the `-a` option which prints the sieve array after sieving is complete. For example, below is a hexidecimal encoding of the donut sieve array for primes in the first octant with norm up to 12000.
 ```shell script
$ ./gintsieve 12000 -a -d
                                                               d34cdbff ffffffec                   
                                                      0aedf3ff c0120c78 da730b2e ffffffef          
                                             6fdcd7ff 1f08ddf0 eed55590 d18b13c7 ffeebdb2          
                                    f65cffef 70113512 369847a3 36968225 e82622d5 1149408b fffffffd 
                           43eedbff f8cba42e 51413f89 082e18b5 88face6d 2a0c2722 5316d330 fffff844 
                  0bffdfef 0e47b6da 25f0b9da f9245714 87cac0ea 8894d858 8f4c51d5 9602fc1a ff54b816 
         f3ceffef ac979b9b be4c61d3 94f6cea5 64332c6d 660ea2b6 1196506d 40fdb31f 2a1e034f f5c0ac56 
5fffffee bb73ba7f c9767527 1739d9f0 c3a99523 bf5e23a4 199f8080 f06afb52 554124ac 4f237261 48e00f02 
```
To see the primality of a 15 x 15 block of Gaussian integers far out the in complex plane, we run:
```shell script
$ ./gintsieve 100000000 300000000 15 15 -a -b
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

Once the python module `gintsieve` is [built](#install), we can import it into python and use it as in the following examples.


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

See [this Jupyter notebook](python/examples.ipynb) for more usage examples.


## Algorithm

With prime generating sieves, we often work with a *sieve array* A and a set of *small primes* P. For each prime p in P, sieving proceeds by crossing off multiples of p in the sieve array A. To implement this on a computer, we define an array of booleans indexed by elements of A. The state of each boolean registers if the underlying index has yet been crossed off by some prime p in P.

In broad strokes, sieving in Z[i] involves the following steps.

1. Initialize a 2-dimensional array A with index (a, b) corresponding to the Gaussian integer a + bi.  In this way, we can view this array as representing a finite subset of the Gaussian integers.  Multiplication by i gives a 4-fold rotational symmetry in the complex plane; complex conjugation (in conjunction with multiplication by i) gives a reflectional symmetry about the line y = x in the complex plane. Such symmetries allow us to reduce the size of the sieving array if it contains redundancies encapsulated by these symmetries.

2. Generate a set of *small primes* P either through calling another sieve or by using the existing sieve array to locate uncrossed elements of small norm. The set P should contain all primes with norm up to sqrt(B), where B is the norm of the largest Gaussian integer in A. In either case, for each prime a + bi in P, leave the sieve array A uncrossed at a + bi but mark off every Gaussian integer multiple of a + bi. Continue this process for each prime in P.

3. The unmarked entries that remain in the array A are exactly the Gaussian primes in A.


[See here](assets/algorithm.pdf) for a discussion of different approaches to generating primes in Z[i] and some motivation for our choice of algorithm.


### Crossing off multiples

Given a Gaussian integer a + bi and a sieve array A (in this implementation, A corresponds to Gaussian integers within a sector or rectangle), we are tasked with crossing off multiples of a + bi from A. Said differently, we must cross off any sieve array element with index corresponding to the complex number (a + bi)(c + di) for c, d in Z. Expanding this, we seek indices of the form (u, v) where
                u = ac - bd,  and  v = ad + bc.
This simply requires a double for-loop (once over c and once over d). Moreover, we can avoid costly multiplications by adding or subtracting c's or d's from the current index u and v while iterating. In the language of abstract algebra, we are crossing off elements of the Z-module generated by (a, b) and (-b, a).

### Donut Sieve

In the rational integers Z, the [wheel sieve](https://en.wikipedia.org/wiki/Wheel_factorization) allows for a more efficient way to manage and step through the sieve array. Primes other than 2, 3, and 5 must lie in one of the 8 distinct residue classes 1, 7, 11, 13, 17, 19, 23, 29 mod 30. In light of this observation, the wheel sieve (with a wheel size of 30) reduces both the memory and complexity of sieving by only tracking these 8 distinguished residue classes. If implemented carefully, a wheel sieve will not add significant overhead to the actual sieving process, and it will reduce to the total size of the sieve array to 8 / 30 = 27% of the original size.

We extend this idea of a wheel sieve to the Gaussian integers. Instead of a wheel corresponding to the quotient Z/30Z rolling over the integers, we consider the finite set of remainders in Z[i]/10Z[i]. Mathematically, this quotient is a **donut** or torus. Just as the primes 2, 3, 5 in the rational integers give rise to a mod-30 wheel, the mod-10 donut arises from the primes 1 + i, 1 + 2i and 2 + i. The first of these three primes sits above the ramifying prime 2, and the second both sit above the splitting prime 5.
- Every Gaussian integer z satisfying z = 0 (mod 1 + i) should be omitted from the sieve array. The prime 1 + i has norm 2, and so one out of every two Gaussian integers will be divisible by 1 + i. In this way, including 1 + i in the donut-base will reduce the sieve array storage and memory by 50%.
- We deal with the primes 1 + 2i and 2 + i together (by the CRT, this is equivalent to working in Z[i]/5Z[i]). Every Gaussian integer satisfying z = 0 (mod 1 + 2i) **or** z = 0 (mod 2 + i) should be omitted from the sieve array. By the inclusion-exclusion principle and the Chinese remainder theorem, there should be exactly 9 solutions to these equations in Z[i]/5Z[i]. The Gaussian integer 5 has norm 25, and so 9 out of every 25 Gaussian integers will be removed from the sieve array. This provides a savings of 9/25 = 36%.
- Combining both of these local conditions, our donut allows us to only consider 32 of the 100 residue classes in Z[i]/10Z[i].

In a wheel sieve, adding more primes to the wheel-base will allow more savings within the sieve array at the expense of increasing the combinatorial complexity of the sieving process itself. Although we stopped with a donut of size 10 for hardware-related reasons (see [implementation](c++-implementation), additional primes to consider as part of a donut-base include:
- using the inert prime 3 in the donut would give a savings of 1/9 = 11% in sieve array size, and
- using the primes 3 + 2i and 3 + 2i, the primes dividing 13, would give a savings of 25/169 = 15% in sieve array size.

In the classes `OctantDonutSieve` and `BlockDonutSieve` we implement the basic structures of this donut sieve. In particular, just as one hard codes gaps between wheel spokes in the wheel sieve, these classes contain variables and methods for accessing and dealing with the donut. As a simple example of how to deal with donut difficulties, the array below shows the 32 distinguished residue classes (with entries corresponding to horizontal gaps) in a mod-10 donut.
```
{{0, 2, 0, 4, 0, 0, 0, 2, 0, 2},
 {4, 0, 0, 0, 2, 0, 4, 0, 0, 0},
 {0, 0, 0, 2, 0, 2, 0, 6, 0, 0},
 {2, 0, 6, 0, 0, 0, 0, 0, 2, 0},
 {0, 4, 0, 0, 0, 4, 0, 0, 0, 2},
 {0, 0, 2, 0, 2, 0, 2, 0, 4, 0},
 {0, 4, 0, 0, 0, 4, 0, 0, 0, 2},
 {2, 0, 6, 0, 0, 0, 0, 0, 2, 0},
 {0, 0, 0, 2, 0, 2, 0, 6, 0, 0},
 {4, 0, 0, 0, 2, 0, 4, 0, 0, 0}}
```

### Segmented Sieve

In a segmented sieve, we break up the sieving array into smaller subsets called **segments**. With these segments in place, we perform sieving (crossing off multiples of small primes) segment by segment. There are advantages and disadvantages to this segmented approach.
- The main disadvantage of segmentation is some additional overhead. Using the following code snippet as an example,
```python
for segment in sieve_array:
    for p in small_primes:
        cross_off_multiples(segment, p)
```
for some primes `p` there will be few or no calls to the function `cross_off_multiples`. As a rule of thumb, each segment should contain at least one multiple of each small prime `'`.
- The main advantage of segmentation is that each segment can fit more readily into computer hardware such as the L2 cache. This will allow the inner for-loop in the snippet above to run more quickly because each segment can be accessed quickly.

In this project, segmentation can be achieved by calling methods in the `BlockSieve` class. In `BlockMoat` and `SegmentedMoat`, we take this approach to explore Gaussian primes.




## C++ Implementation

The aforementioned algorithm is implemented in a C++ library. `BaseSieve` is an abstract base class with some basic sieving methods. Classes derived from this include `OctantSieve`, `OctantDonutSieve`, `SectorSieve`, `BlockSieve`, and `BlockDonutSieve`. Each derived class has its own method for initiating and accessing the sieve array. See the [usage examples](#command-line-usage) for various text representations of these sieve arrays.


Results of this repository can be verified with [http://oeis.org/A091100](http://oeis.org/A091100) to check for accuracy. These tests can be performed by running the following.
```shell script
$./ginttest
```

In classes `OctantSieve`, `SectorSieve`, and `BlockSieve`, the 2-dimensional sieve array A is stored with a `vector<vector<bool>>` container. This C++ object is memory-efficient: each boolean is stored as a single bit in memory (this in itself gives a 224-fold savings over python 3.7 which requires 28 bytes of memory to store a boolean value).

In implementing the donut sieve in the classes `OctantDonutSieve` and `BlockDonutSieve`, each 10 x 10 block of Gaussian integers corresponds to a full *donut roll*. The [donut sieve](#donut-sieve) requires holding 32 residue classes for every 10 x 10 block of Gaussian integers. Said differently, every 10 x 10 block of Gaussian integers requires 32 bits of information to store its current state in the sieve process. Conveniently, a C++ `int` typically also requires 32 bits of memory space. In this way, in donut-based classes , our sieve array is a `vector<vector<int>>` container in C++.



## Applications

This library can be used to generated data pertaining to several unsolved problems in number theory.
- [Gaussian prime races in sectors](#gaussian-prime-races-in-sectors)
- [Angular distribution of Gaussian primes](#angular-distribution-of-gaussian-primes)
- [The Gaussian moat problem](#the-gaussian-moat-problem)

All three of these applications have implementations within the [Python bindings](#python-bindings) and can be found inside the `gintsieve` python module. [This Jupyter notebook](python/examples.ipynb) gives some examples.


### Gaussian prime races in sectors
Prime number races have been thoroughly studied beginning with observation of Chebyshev. Odd prime numbers fall into two disjoint camps: primes of the form 4k + 1, and primes of the form 4k + 3. Chebyshev noted that when counting the number of primes in these two camps up to some threshold x, the primes of the form 4k + 3 often seemed to lead the race. This phenomenon has been extensively quantified and generalized in the realm of the rational integers, and is known as *Chebyshev bias*.

In the Gaussian integers, one analog of 4k + 1 and 4k + 3 primes is considering Gaussian primes by sector. A [classical result](http://gdz.sub.uni-goettingen.de/dms/resolveppn/?PPN=GDZPPN002365162) of Hecke is that the angles determined by Gaussian primes are equidistributed. In other words, given two sectors in the complex plane with equal central angles, Hecke proved that both sectors have the same asymptotic number of primes as the radius of those two sectors grows large. Just as Chebyshev raced the primes in the residue classes 4k + 1 and 4k + 3, one can hold a *fair race* among Gaussian primes in two sectors with equal central angle.

The `SectorSieve` class performs sieving in a specified sector in the complex plane. This can be used to generate data and observe Chebyshev biases in Gaussian prime races.


### Angular distribution of Gaussian primes
Instead of restricting to two-way prime number races, one can consider races with many participants. As a first example of this in the rational primes, odd primes fall into one of the four camps 8k + 1, 8k + 3, 8k + 5, and 8k + 7. As a consequence, one can consider the four-way race between these four disjoint residue classes mod 8. This many-way race can just as easily be considered with Gaussian primes. In particular, a collection of sectors having the same central angle, one can consider the distribution of the counts of Gaussian primes within these sectors as the bound on the sectors' radii grows.


### The Gaussian moat problem
In the Gaussian moat problem, imagine a frog jumping between lily pads situated in the complex plane at exactly the Gaussian primes. The frog, having finite strength, can only take jumps of bounded size. A famous unsolved problems asks if such a frog could jump arbitrarily far from its starting point. If this frog fails to do so, it must have encountered some *Gaussian moat* preventing it from making progress.

This problem, lying at the interface between number theory and percolation theory, has been extensively studied since it was first posed in 1962. In 2004, Tsuchimura computed the state of the result by showing that a frog starting at the origin and taking jumps of distance at most 6 will encounter a Gaussian moat. It is a folklore conjecture that given any bounded jump size, a moat will always be encountered.

Beginning with a bound on the frog's jump, form a graph whose vertices the Gaussian primes. Two vertices are adjacent if the distance between them is less than or equal to the jump bound. Three classes within this library explore the *main connected component* of this graph, that is, the connected component containing the Gaussian prime 1 + i (the prime closest to the origin).


- In the class `OctantMoat`, the main connected component is explored using a depth-first search approach. Initially, all primes are generated and held in memory using the `OctantSieve` class. Once the sieve array is larger than the memory capabilities of the computer calling this class, this approach fails.
- In `SegmentedMoat`, the same component is explored using a segmented sieving approach. Instances of `BlockSieve` are called and explored; boundary data is passed from one instance to the next. This approach is to similar to the algorithm in Tsuchimura's paper *Computational Results for Gaussian Moat Problem*. Only the total component size is returned by this algorithm.
- In `BlockMoat`, we consider a vertical strip in the complex plane running between the real-axis and the diagonal line defining the top of the first octant. As in `SegmentedMoat`, instances of `BlockMoat` are called and boundary data is shared. This algorithm searches for a sequence of adjacent Gaussian primes spanning the vertical strip. If no such sequence is found, then a moat is present, and as a consequence, the main component is finite. A similar approach was taken in the paper *A stroll Through the Gaussian Primes* by Gethner, Wagon, and Wick.

All of these moat-exploring classes can be accessed from the `gintmoat` executable. Command line options include
```
Usage: ./gintmoat jumpSize [realPart] [option1] [option2] ...
Explore the connected component at the origin in the Gaussian moat problem.
    jumpSize            The jump bound giving adjacency relation among primes.
    realPart            The real-part of the vertical strip to explore if in
                        vertical mode.

Exploration modes:
    --origin            Explore the connected component starting at the origin until
                        an impassable moat is encountered. Default exploration mode.
    --segmented         Use a segmented approach to explore the connected component.
                        Algorithm is similar to that in Tsuchimura paper. This approach
                        only counts the size of the component.
    --vertical          Search for a Gaussian moat along a thin vertical strip starting
                        at real-part x. Used to show a component is finite.

Options:
    -h, --help          Print this help message.
    -v, --verbose       Display progress.
    -p, --printprimes   Print the real and imag part of primes in the connected component
                        if in origin mode.
```

For example, to print all primes that can be reach with jumps up to distance 1.5, we run:
```shell script
$ ./gintmoat 1.5 -p

1 1
2 1
3 0
3 2
4 1
5 2
6 1
7 0
7 2
8 3
9 4
8 5
10 3
11 4
The main connected component has size: 14
The farthest out prime in component has coordinates: 11 4
```
To get the size of the component with jumps at most distance 4.5, we could run the following (and wait nearly an hour).
```shell script
$ ./gintmoat 4.5 --segmented


The main connected component has size: 273791623
```

Many other graph theoretic statistics could be gathered using a similar implementation. Examples of interesting data to consider include:
- The distribution of the vertex degrees in the main component
- Number and location of bridges
- Size of other connected components (away from the origin)
- "Lakes" and other regions void of primes

This project could readily be adapted to explore such phenomena.


## License

This project is released under the [MIT license](https://opensource.org/licenses/MIT).