# cython: language_level=3
# distutils: language = c++

from gaussian_prime_sieve cimport cythonTest2

def p2(x):
    return cythonTest2(x)

def hello():
    print("Hello World")

def foo():
    print('bar')