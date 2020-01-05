# cython: language_level=3
# distutils: language = c++

from gaussian_prime_sieve cimport count

def get_primes_to_norm(x):
    return count(x)

def hello():
    print("Hello World")