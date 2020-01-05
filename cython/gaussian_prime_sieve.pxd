# Read this whole section
# https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html#add-public-attributes
cdef extern from '../cpp/Main.cpp':
    pass


cdef extern from '../cpp/Main.hpp':
    int cythonTest2(int x)




