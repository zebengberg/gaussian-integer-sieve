# cython: language_level=3



from gintsieve_externs cimport gPrimes, gPrimesCount



# cython casts vector to a list
# see https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html#standard-library


cpdef gprimes(long x):
    return gPrimes(x)


cpdef gprimes_count(long x):
    return gPrimesCount(x)