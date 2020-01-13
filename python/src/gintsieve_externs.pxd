from libcpp.vector cimport vector
from libcpp.pair cimport pair
cimport numpy as np


cdef extern from '../../include/cython_bindings.hpp':
    vector[pair[long, long]] gPrimes(long)
    vector[pair[long, long]] gPrimesSegment(long, long, long)
    unsigned long gPrimesCount(long)
    unsigned long gPrimesSegmentCount(long, long, long)
    pair[long long *, unsigned long] gPrimesArray(long);



