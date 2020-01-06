from libcpp.vector cimport vector
from libcpp.pair cimport pair


cdef extern from '../cpp/main.hpp':
    vector[pair[long, long]] gPrimes(long)
    unsigned long gPrimesCount(long)
    vector[pair[long, long]] gPrimesRegion(long, long, long)


