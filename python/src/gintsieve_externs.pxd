from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t, uint64_t
cimport numpy as np


cdef extern from '../../include/cython_bindings.hpp':
    vector[pair[uint32_t, uint32_t]] gPrimesToNorm(uint64_t)
    vector[pair[uint32_t, uint32_t]] gPrimesInSector(uint64_t, double, double)
    vector[pair[uint32_t, uint32_t]] gPrimesInBlock(uint32_t, uint32_t, uint32_t, uint32_t)

    uint64_t gPrimesToNormCount(uint64_t)
    uint64_t gPrimesInSectorCount(uint64_t, double, double)
    uint64_t gPrimesInBlockCount(uint32_t, uint32_t, uint32_t, uint32_t)

    # Couldn't figure out how to use a pointer to uint32_t.
    # Instead using unsigned int, which is a 32 bit integer in cython
    pair[unsigned int *, uint64_t] gPrimesToNormAsArray(uint64_t)
    pair[unsigned int *, uint64_t] gPrimesInSectorAsArray(uint64_t, double, double)
    pair[unsigned int *, uint64_t] gPrimesInBlockAsArray(uint32_t, uint32_t, uint32_t, uint32_t)
