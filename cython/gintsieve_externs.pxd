# Read this whole section
# https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html#add-public-attributes

# cdef extern from '../cpp/BaseSieve.cpp':
#     pass

# cdef extern from '../cpp/QuadrantSieve.cpp':
#     pass

#cdef extern from '../cpp/main.cpp':
#    pass





cdef extern from '../cpp/main.hpp':
    long count(long x)




