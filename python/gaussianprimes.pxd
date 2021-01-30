from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t, uint64_t, int32_t
cimport numpy as np

# work around for bug with pointers
# see https://github.com/cython/cython/issues/534
ctypedef int32_t * intptr


cdef extern from 'cython_bindings.hpp':
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

  # Functions accessing moat data
  pair[intptr, uint64_t] moatMainComponent(double)
  vector[pair[intptr, uint64_t]] moatComponentsToNorm(double, uint64_t)
  vector[pair[intptr, uint64_t]] moatComponentsInBlock(double, int32_t, int32_t, int32_t, int32_t)