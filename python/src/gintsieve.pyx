# cython: language_level=3

import numpy as np
cimport numpy as np
from cython cimport view
import matplotlib.pyplot as plt
from libc.stdint cimport uint32_t, uint64_t

from gintsieve_externs cimport gPrimesToNorm, gPrimesInSector, gPrimesInBlock,\
    gPrimesToNormCount, gPrimesInSectorCount, gPrimesInBlockCount,\
    gPrimesToNormAsArray, gPrimesInSectorAsArray, gPrimesInBlockAsArray,\
    angularDistribution

cpdef gprimes(uint64_t x):
    """Generate a GintList object of primes up to norm x."""
    # Cython casts vector to an array.
    return GintList(gPrimesToNorm(x), x)

cpdef gprimes_sector(uint64_t x, double alpha, double beta):
    """Generate a GintList object of primes up to norm x in given sector."""
    return GintList(gPrimesInSector(x, alpha, beta), x, alpha, beta)

cpdef gprimes_block(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy):
    """Generate a GintList object of primes up to norm x in given sector."""
    return GintList(gPrimesInBlock(x, y, dx, dy), x, y, dx, dy)

cpdef count_gprimes(uint64_t x):
    """Count primes up to norm x."""
    return gPrimesToNormCount(x)

cpdef count_gprimes_sector(uint64_t x, double alpha, double beta):
    """Count primes up to norm x."""
    return gPrimesInSectorCount(x, alpha, beta)

cpdef count_gprimes_block(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy):
    """Count primes up to norm x."""
    return gPrimesInBlockCount(x, y, dx, dy)

cpdef gprimes_as_np(uint64_t x):
    both = gPrimesToNormAsArray(x)
    cdef uint32_t *ptr = both.first
    cdef uint64_t size = both.second
    # Casting c++ pointer returned by gPrimesToNormArray to a memory view object
    cdef view.array primes = <np.uint32_t[:size]> ptr
    cdef np_primes = np.asarray(primes)
    return np_primes

cpdef gprimes_sector_as_np(uint64_t x, double alpha, double beta):
    both = gPrimesInSectorAsArray(x, alpha, beta)
    cdef uint32_t *ptr = both.first
    cdef uint64_t size = both.second
    # Casting c++ pointer returned by gPrimesToNormArray to a memory view object
    cdef view.array primes = <np.uint32_t[:size]> ptr
    cdef np_primes = np.asarray(primes)
    return np_primes

cpdef gprimes_block_as_np(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy):
    both = gPrimesInBlockAsArray(x, y, dx, dy)
    cdef uint32_t *ptr = both.first
    cdef uint64_t size = both.second
    # Casting c++ pointer returned by gPrimesToNormArray to a memory view object
    cdef view.array primes = <np.uint32_t[:size]> ptr
    cdef np_primes = np.asarray(primes)
    return np_primes




class GintList(list):
    """Class to wrap sieved lists."""
    def __init__(self, gints, x, x1=-1, x2=-1, x3=-1):
        list.__init__(self, gints)
        self.x = x
        if x3 != -1:
            self.y = x1
            self.dx = x2
            self.dy = x3
            self.sieve = 'block'
        elif x2 != -1:
            self.alpha = min(x1, x2)
            self.beta = max(x1, x2)
            self.sieve = 'sector'
        else:
            self.sieve = 'to_norm'

    def visualize(self, full_disk=False, save=False):
        """Plot Gaussian primes with Matplotlib."""

        reals = [p[0] for p in self]
        imags = [p[1] for p in self]

        plt.subplots(figsize=(8, 8))

        if self.sieve == 'to_norm':
            # Drawing x and y axes.
            plt.axhline(0, color='red')
            plt.axvline(0, color='red')

            if full_disk:
                # Filling out second quadrant.
                reals, imags = reals + [-t for t in imags], imags + reals

                # Filling out third and fourth quadrants.
                reals += [-t for t in reals]
                imags += [-t for t in imags]

                plt.plot(reals, imags, 'bo', markersize=100 / (self.x ** .5))

            else:
                plt.plot(reals, imags, 'bo', markersize=200 / (self.x ** .5))
        else:
            # Plotting segment.
            plt.plot(reals, imags, 'ro', markersize=300 / self.z)

        if save:
            plt.savefig('sieve_visual.png')

        plt.show()


    def to_complex(self):
        """Convert list of tuples to list of complex numbers."""
        return [complex(pair[0], pair[1]) for pair in self]

    def to_np(self):
        """Convert list of tuples to np.array."""
        return np.array(self).transpose()

    def sector_race(self, a, b, c, d):
        """Gaussian prime race in sectors."""
        if b < a or d < c or a < 0 or c < 0 or b > np.pi/2 or d > np.pi/2:
            raise ValueError('Check your intervals.')
        p = self.to_np()
        norms = p[0] ** 2 + p[1] ** 2
        angles = np.arctan2(p[0], p[1])
        runner1 = ((angles > a) & (angles < b)).cumsum()
        runner2 = ((angles > c) & (angles < d)).cumsum()

        plt.subplots(figsize=(8, 8))
        plt.plot(norms, runner1 - runner2, 'b-')
        plt.title('Gaussian prime race in sectors')
        plt.xlabel('norm')
        plt.ylabel('$\pi({}, {}) - \pi({}, {})$'.format(a, b, c, d))
        plt.axhline(0, color='red')
        plt.show()

cpdef race(a, b, c, d):
    # histogram?
    # plot norm against cum sum?
    # weight them in different ways?
    # specify norm or let function specify it?
    pass

cpdef moat():
    # get prmes in block
    pass

cpdef angular_dist(x, n_sectors):
    """Make histogram of number of Gaussian primes up to norm x in equispaced sectors."""
    data = np.array(angularDistribution(x, n_sectors))
    min = np.percentile(data, 1)
    max = np.percentile(data, 99)
    plt.hist(data, bins=30, range=(min, max))
    plt.show()
    return data




