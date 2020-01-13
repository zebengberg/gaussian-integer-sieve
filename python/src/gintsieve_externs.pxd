from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t, uint64_t
cimport numpy as np


cdef extern from '../../include/cython_bindings.hpp':
    vector[pair[uint32_t, uint32_t]] gPrimes(uint64_t)
    vector[pair[uint32_t, uint32_t]] gPrimesSegment(uint32_t, uint32_t, uint32_t)
    uint64_t gPrimesCount(uint64_t)
    uint64_t gPrimesSegmentCount(uint32_t, uint32_t, uint32_t z)
    pair[unsigned int *, uint64_t] gPrimesAsArray(uint64_t)