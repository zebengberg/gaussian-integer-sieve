# gaussian-prime-sieve

This is an implementation of sieve of Eratosthenes to generate primes in the Gaussian integers.  The Gaussian integers $\mathbb{Z}[i]$ are complex numbers of the form $a+bi$ for integers $a$ and $b$, and enjoy many of the same properties of the rational integers $\mathbb{Z}$.  In particular, Gaussian integers factorize uniquely into primes.  Similar to the situation in $\mathbb{Z}$, sieving algorithms can be used to compute primes in $\mathbb{Z}[i]$.

The goal of this project is to build an efficient implementation of a prime-generating sieve in the Gaussian integers in both Python and C, and to provide a library for deeper exploration.


## Examples

```Python
>>>  from simple_sieve import simple_sieve, visualize

# Generate Gaussian primes with norm up to 50 sorted lexicographically.
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

# Generate Gaussian primes with norm up to 50 sorted by norm.
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
![First Quadrant](https://github.com/zebengberg/gaussian-prime-sieve/blob/master/images/first_quadrant.png)

```python
# Plot the Gaussian primes in the full complex plane with norm up to 1000.
>>>  x = 1000
>>>  P = simple_sieve(x)
>>>  visualize(P, x, full=True)
```
![](full_plane.png)



## Algorithm

Below is a sketch of this implementation of the sieve of Eratosthenes in the Gaussian integers.  Suppose we want to generate Gaussian primes with norm up to $x$, some large real number.

1. Initialize an empty 2-dimensional array with index $(a, b)$ corresponding to the Gaussian integer $a + bi$.  In this way, we can simply view this array as a finite subset of the Gaussian integers.  By restricting to the first quadrant of the complex plane, we take $a$ and $b$ to be nonnegative integers with $a^2 + b^2 \le x$.  Mark off the array at the Gaussian integers $0, 1$, and $i$.  These Gaussian integers are neither prime nor composite.
2. Locate the unmarked Gaussian integer $\pi$ with smallest norm.  Necessarily $\pi$ is prime.  Leave the array empty at $\pi$ but mark off every Gaussian integer multiple of $\pi$.
3. Loop over the process in step 2 by marking off multiples of the other previously unmarked Gaussian primes.  Continue looping until every unmarked entry with norm less than $\sqrt x$ has been used for sieving.
4. The unmarked entries that remain in the array are exactly the Gaussian primes.

Compared with the situation in $\mathbb{Z}$, step 2 is somewhat more involved in $\mathbb{Z}[i]$.  There are several cases to consider depending on the splitting behavior of the prime $\pi$ in question.


- If $\pi = a + bi$ with $b=0$, necessarily $a = p$ is prime with $p\equiv 3 \bmod 4$.  In this case, $\pi$ is a rational integer, and we can simply mark off all Gaussian integers of the form $(s + it) p$ to eliminate multiples of $\pi$.  Here, the multiples $(s + it) p$ for $s,t\in\mathbb{Z}$ form a square lattice perpendicular to the coordinate axes, and so the marking can be implemented with two nested for-loops.


- If $\pi = a + bi$ with $ab \ne 0$, then the norm $N(\pi) = a^2 + b^2$ must be a prime.  Let $p$ be this prime.  Either $p=2$ or $p\equiv 1\bmod 4$, and $p$ is said to ramify or split, respectively.  In either case, $p$ is a prime of degree 1, meaning that the additive quotient group $\mathbb{Z}[i] / \pi \mathbb{Z}[i]$ is a group of order $p^1$.  Furthermore, the element $\pi$ has order $p$ in the group $\mathbb{Z}[i] / p\mathbb{Z}[i]$.  To construct multiples of $\pi$ in $\mathbb{Z}[i]$, we will construct elements of additive group generated by $\pi$, then translate them by multiplies of $p$.  Specifically, if we let $s + it$ be one of the $p$ distinct elements in the additive subgroup generated by $\pi \bmod p$, and we mark off all Gaussian integers of the form $(s + it)p$.  As in the first case, this marking-off can be implemented with nested for-loops.  The elements $s + it$ can be obtained by reducing the coefficients of $k(a+bi)$ modulo $p$ for $k=0, 1, 2, \dots, p-1$.

The theorem that makes this procedure succeed can be stated as follows.  With $\pi$ a degree 1 prime, given any Gaussian integer $\alpha$, there exists a Gaussian integer $\delta$ and a rational integer $d$ such that $\pi\alpha = d\pi + \delta p.$  In other words, the Gaussian-integer-multiple of $\pi$ given by $\alpha \pi$ can be obtained through this aforementioned process.

There are several ways to prove this.  One explicit way to do so is to let $\pi = a + bi$ and $\alpha = u + vi$.  Let $c \in \mathbb{Z}$ be the multiplicative inverse of $b$ modulo $p$.  Then $bc \equiv 1\bmod p$, and so $bc \equiv 1\bmod \overline{\pi}$ since $\overline{\pi}$ divides $p$.  Because $a\equiv bi\bmod \overline{\pi}$, it follows that $ac \equiv i\bmod \overline{\pi}$.  Therefore
$$u + vi \equiv u + vac \bmod \overline{\pi},$$
which shows that $\alpha \equiv d\bmod \overline{\pi}$ for some rational integer $d$.  Unpacking this congruence, we have
$$\alpha  -d = \overline{\pi} \delta$$
for some $\delta \in\mathbb{Z}[i]$.  Multiplying both sides of this equation by $\pi$ gives the desired conclusion.








## License

This project is licensed under the terms of the MIT license.
