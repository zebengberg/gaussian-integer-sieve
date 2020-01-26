from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t, uint64_t, int32_t
cimport numpy as np

# Work around for bug; see https://github.com/cython/cython/issues/534
ctypedef uint32_t * intptr


cdef extern from '../../include/cython_bindings.hpp':
    vector[pair[uint32_t, uint32_t]] gPrimesToNorm(uint64_t)
    vector[pair[uint32_t, uint32_t]] gPrimesInSector(uint64_t, double, double)
    vector[pair[uint32_t, uint32_t]] gPrimesInBlock(uint32_t, uint32_t, uint32_t, uint32_t)

    uint64_t gPrimesToNormCount(uint64_t)
    uint64_t gPrimesInSectorCount(uint64_t, double, double)
    uint64_t gPrimesInBlockCount(uint32_t, uint32_t, uint32_t, uint32_t)

    # Couldn't figure out how to use a pointer to uint32_t inside pair
    # Instead using unsigned int, which is a 32 bit integer in cython
    pair[intptr, uint64_t] gPrimesToNormAsArray(uint64_t)
    pair[intptr, uint64_t] gPrimesInSectorAsArray(uint64_t, long double, long double)
    pair[intptr, uint64_t] gPrimesInBlockAsArray(uint32_t, uint32_t, uint32_t, uint32_t)

    vector[uint64_t] angularDistribution(uint64_t, uint32_t)

    # Using this class to transfer race data to numpy
    cdef cppclass SectorRace:
        SectorRace() except +
        SectorRace(uint64_t, uint64_t, long double, long double, long double, long double) except +
        pair[intptr, uint64_t] getFirstSector()
        pair[intptr, uint64_t] getSecondSector()
        int32_t * getNormData()
