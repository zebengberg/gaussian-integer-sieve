from libcpp.vector cimport vector
from libcpp.pair cimport pair


cdef extern from '../cpp/main.hpp':
    vector[pair[long, long]] gPrimes(long)
    vector[pair[long, long]] gPrimesSegment(long, long, long)
    unsigned long gPrimesCount(long)



