from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t, uint64_t
cimport numpy as np


cdef extern from '../../include/cython_bindings.hpp':
    vector[pair[uint32_t, uint32_t]] gPrimes(uint64_t)
    vector[pair[uint32_t, uint32_t]] gPrimesSegment(uint32_t, uint32_t, uint32_t)
    uint64_t gPrimesCount(uint64_t)
    uint64_t gPrimesSegmentCount(uint32_t, uint32_t, uint32_t z)
    # Couldn't figure out how to use a pointer to uint32_t.
    # Instead using unsigned int, which is a 32 bit integer in cython
    pair[unsigned int *, uint64_t] gPrimesAsArray(uint64_t)