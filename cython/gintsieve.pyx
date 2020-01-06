# cython: language_level=3

from gintsieve_externs cimport count

def c(x):
    return count(x)

def hello():
    print("Hello World")