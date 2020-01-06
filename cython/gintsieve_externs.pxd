from libcpp.vector cimport vector
from libcpp.pair cimport pair


cdef extern from '../cpp/main.hpp':
    vector[pair[long, long]] gPrimes(long x)
    long gPrimesCount(long x)
    vector[pair[long, long]] gPrimesRegion(long x, long y, long z)


