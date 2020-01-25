# cython: language_level=3

import numpy as np
cimport numpy as cnp
from cython cimport view
import matplotlib.pyplot as plt
from libc.stdint cimport uint32_t, uint64_t, int32_t

from gintsieve_externs cimport gPrimesToNorm, gPrimesInSector, gPrimesInBlock,\
    gPrimesToNormCount, gPrimesInSectorCount, gPrimesInBlockCount,\
    gPrimesToNormAsArray, gPrimesInSectorAsArray, gPrimesInBlockAsArray,\
    angularDistribution, sectorRace

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
    cdef view.array primes = <cnp.uint32_t[:size]> ptr
    cdef np_primes = np.asarray(primes)
    # De-flattening array.
    np_primes = np_primes.reshape(size // 2, 2)
    return np_primes.transpose()

cpdef gprimes_sector_as_np(uint64_t x, double alpha, double beta):
    both = gPrimesInSectorAsArray(x, alpha, beta)
    cdef uint32_t *ptr = both.first
    cdef uint64_t size = both.second
    # Casting c++ pointer returned by gPrimesInSectorAsArray to a memory view object
    cdef view.array primes = <cnp.uint32_t[:size]> ptr
    cdef np_primes = np.asarray(primes)
    # De-flattening array.
    np_primes = np_primes.reshape(size // 2, 2)
    return np_primes.transpose()

cpdef gprimes_block_as_np(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy):
    both = gPrimesInBlockAsArray(x, y, dx, dy)
    cdef uint32_t *ptr = both.first
    cdef uint64_t size = both.second
    # Casting c++ pointer returned by gPrimesInBlockAsArray to a memory view object
    cdef view.array primes = <cnp.uint32_t[:size]> ptr
    cdef np_primes = np.asarray(primes)
    # De-flattening array.
    np_primes = np_primes.reshape(size // 2, 2)
    return np_primes.transpose()




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


cpdef angular_dist(uint64_t x, uint32_t n_sectors):
    """Make histogram of number of Gaussian primes up to norm x in equispaced sectors."""
    data = np.array(angularDistribution(x, n_sectors))
    min = np.percentile(data, 1)
    max = np.percentile(data, 99)
    plt.hist(data, bins=30, range=(min, max))
    plt.show()
    return data

class SectorRace:
    """Class to hold data on Gaussian prime races."""
    def __init__(self, x, a, b, c, d, n_bins=1000):
        if (b - a) - (d - c) > 0.0000001:
            raise ValueError('Unfair race; try again with two equal-sized intervals.')
        for angle in [a, b, c, d]:
            if angle < 0 or angle > np.pi / 4:
                raise ValueError('Stay inside the first octant.')

        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.x = x

        self.n_bins = n_bins
        self.bin_width = self.x // self.n_bins
        self.norms = np.linspace(self.bin_width, self.x, self.n_bins)
        normalizer = lambda x: np.sqrt(x) / np.log(x)
        self.normalize = normalizer(self.norms)



        cdef int32_t *ptr = sectorRace(x, a, b, c, d, n_bins)
        cdef view.array norm_data = <cnp.int32_t[:n_bins]> ptr
        self.norm_data = np.asarray(norm_data)




    def plot_race(self, normalize=True):
        """Plot """
        plt.subplots(figsize=(8, 8))
        plt.plot(self.norms, self.norm_data, 'b-')
        if normalize:
            plt.plot(self.norms, self.normalize, 'r-')
            plt.plot(self.norms, -self.normalize, 'r-')

        plt.title('Gaussian prime race in sectors')
        plt.xlabel('$x$')
        plt.ylabel('$\pi(x; {}, {}) - \pi(x; {}, {})$'.format(self.a, self.b, self.c, self.d))
        plt.axhline(0, color='red')
        plt.show()

    def hist(self, bins='auto'):
        """Plot Shanks-style histogram of race progress."""
        shanks = self.norm_data / self.normalize
        plt.hist(shanks, bins=bins)
        plt.show()

    def density(self):
        """Approximate log-density of how often first sector leads on interval [1, x]."""
        leader = np.where(self.norm_data > 0, 1, 0)
        ties = np.where(self.norm_data == 0, 0.5, 0)
        weights = np.arange(1, self.n_bins + 1)
        weights = 1 / weights
        d = np.sum(weights * (leader + ties))
        d /= np.sum(weights)
        return d

