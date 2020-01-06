# cython: language_level=3

from gaussian_prime_sieve cimport count

def c(x):
    return count(x)

def hello():
    print("Hello World")