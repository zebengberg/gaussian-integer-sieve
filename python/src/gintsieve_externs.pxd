from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t, uint64_t, int32_t
cimport numpy as np

# Work around for bug; see https://github.com/cython/cython/issues/534
ctypedef int32_t * intptr


cdef extern from '../../include/cython_bindings.hpp':
    # TODO: possibly delete these three if don't want to offer native python support
    vector[pair[int32_t, int32_t]] gPrimesToNorm(uint64_t)
    vector[pair[int32_t, int32_t]] gPrimesInSector(uint64_t, double, double)
    vector[pair[int32_t, int32_t]] gPrimesInBlock(uint32_t, uint32_t, uint32_t, uint32_t)

    uint64_t gPrimesToNormCount(uint64_t)
    uint64_t gPrimesInSectorCount(uint64_t, double, double)
    uint64_t gPrimesInBlockCount(uint32_t, uint32_t, uint32_t, uint32_t)

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
        intptr getNormData()

    # Using this class to transfer moat data to numpy
    cdef cppclass OctantMoat:
        OctantMoat() except +
        OctantMoat(uint64_t, double) except +
        void exploreComponent(int32_t, int32_t)
        void exploreAllComponents()
        pair[intptr, uint64_t] getCurrentComponent()
        pair[intptr, uint64_t] getUnexplored()
        vector[pair[intptr, uint64_t]] getAllComponents()
