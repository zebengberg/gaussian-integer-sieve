# gaussian-prime-sieve

>A program to generate primes in the Gaussian integers using the Sieve of Eratosthenes.


## Table of Contents

- [Background](#background)
- [Install](#install)
- [Examples](#examples)
- [Algorithm](#algorithm)
- [License](#license)


## Background

The **Gaussian integers** are complex numbers of the form a + bi for integers a and b. Here i is the **imaginary unit** with the property that i * i = -1. The set of Gaussian integers, denoted by Z[i], forms a mathematical structure that enjoys many of the same properties as those of the rational integers, denoted by Z. Gaussian integers can be factored and the notion of primality is well-defined. The **norm** of a Gaussian integer gives a measure of its size; the analogy in the rational integers is absolute value. With these tools in hand, sieving ideas can be used to compute primes in this extension of our usual integer system.

This project is an implementation of the *Sieve of Eratosthenes* to generate primes in the Gaussian integers. While the Gaussian integers present some difficulties that mostly nonexistent in the context of the rational integers, nearly all of the basic computational sieve ideas can be extended to Z[i]. In particular, both the wheel sieve and the segmented sieve can be extended to the Gaussian integers.

The goal of this project is to build an efficient implementation of a prime generating sieve in the Gaussian integers in both Python and C++, and to provide a library for deeper exploration. The [*primesieve*](https://github.com/kimwalisch/primesieve) project implements the current state of the art prime generator for rational integers. *primesieve* includes both algorithmic optimizations as well as hardware considerations needed in order to improve runtime performance. In this repository, I hope to extend some of the more accessible algorithmic and hardware-optimization ideas from *primesieve*.


## Install

This entire repository can be cloned or downloaded for use.

The python module `simple_sieve.py` can be imported into a python console or run directly from the command line. After an initial list of "small" primes is generated and stored to file, the C++ file `simple_sieve.cpp` can be compiled and run. See the first example below for details.

## Examples

To generate all Gaussian integer primes with norm up to 10000, enter
```shell script
$ python simple_sieve.py 10000
```
into the command line from within the `python` folder of this repository. This command will write the primes generated to a text file which can then be read by the C++ program `simple_sieve.cpp`. A text file containing a list of primes with norm up to sqrt(x) is needed in order to use `simple_sieve.cpp` to generate primes with norm up to x. 

The python module `simple_sieve.py` can also be imported into a python script as seen in the examples below.

```Python
>>>  from simple_sieve import simple_sieve, visualize

# Generate Gaussian integer primes with norm up to 50 sorted lexicographically.
>>>  simple_sieve(50)
[(1, 1),
 (1, 2),
 (1, 4),
 (1, 6),
 (2, 1),
 (2, 3),
 (2, 5),
 (3, 0),
 (3, 2),
 (4, 1),
 (4, 5),
 (5, 2),
 (5, 4),
 (6, 1),
 (7, 0)]

# Generate Gaussian integer primes with norm up to 50 sorted by norm.
>>>  simple_sieve(50, sort=True)
[(1, 1),
 (1, 2),
 (2, 1),
 (3, 0),
 (2, 3),
 (3, 2),
 (1, 4),
 (4, 1),
 (2, 5),
 (5, 2),
 (1, 6),
 (6, 1),
 (4, 5),
 (5, 4),
 (7, 0)]

# Plot the Gaussian primes in the first quadrant with norm up to 10000.
>>>  x = 10000
>>>  P = simple_sieve(x)
>>>  visualize(P, x)
```
![First Quadrant](/images/first_quadrant.png)

```python
# Plot the Gaussian primes in the full complex plane with norm up to 1000.
>>>  x = 1000
>>>  P = simple_sieve(x)
>>>  visualize(P, x, full_disk=True)
```
![Full Plane](/images/full_plane.png)

## Implementation

With prime generating sieves, we often work with a *sieve array* A and a set of primes P. For each prime p in P, the sieve proceeds by crossing off multiples of p in the sieve array A. To implement this on a computer, we define an array of booleans indexed by elements of A. The state of each boolean tracks if the element has yet been crossed off by some prime p in P.

Working within the Gaussian integers, our sieve array A is a 2-dimensional object indexed by Gaussian integers. In native python, this can be implemented with a list of lists. We can also use a 2-dimensional *numpy* array. In C++, we can use the built-in array, or vector, or another class of container.

In python 3.7, the boolean `True` takes up 28 bytes of memory storage. The vast majority of this is overhead; in theory, a boolean value should only occupy a single bit of space. Even in C++, a single boolean value requires a full byte of memory. One issue involved with storing a boolean value into a smaller storage space is that modern computers do not allow pointers to individual bits in memory. In C++, there are various ways around this issue such as the `bitset` object.

A different way to resolve this issue is to include some of the ideas of the wheel sieve. Primes other than 2, 3, and 5 must lie in one of the 8 distinct residue classes 1, 7, 11, 13, 17, 19, 23, 29 mod 30. In light of this observation, the wheel sieve reduces both the memory and complexity of sieving by only tracking these 8 distinguished residue classes. If implemented carefully, a wheel sieve should not add significant overhead to the actual sieving process. Moreover, this approach dovetails nicely with modern computer architecture: a single byte can track the state of 30 consecutive integers from an integer sieve array, allowing for a 30-fold reduction in memory usage (in a C++ implementation). Extending this wheel sieve to the Gaussian integers, both the memory and complexity of the sieve can be greatly reduced.

In `simple_sieve.py` and `simple_sieve.cpp`, I use native booleans to track the state of each Gaussian integer in the sieve array. This will create significant overhead in both python (up to 28 * 8 = 224-fold waste) and C++ (8-fold).  In our donut sieve (TODO: implement this), I use the aforementioned ideas to reduce and align the storage requirements of the sieve to modern hardware architecture.

Results of this repository can be verified with [http://oeis.org/A091100](http://oeis.org/A091100) to check for accuracy.


## Simple sieving algorithm

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



## License

This project is released under the [MIT license](https://opensource.org/licenses/MIT).