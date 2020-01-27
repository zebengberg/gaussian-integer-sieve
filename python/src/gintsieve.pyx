# cython: language_level=3

import numpy as np
cimport numpy as cnp
from cython cimport view
import matplotlib.pyplot as plt
from libc.stdint cimport uint32_t, uint64_t, int32_t

from gintsieve_externs cimport gPrimesToNorm, gPrimesInSector, gPrimesInBlock,\
    gPrimesToNormCount, gPrimesInSectorCount, gPrimesInBlockCount,\
    gPrimesToNormAsArray, gPrimesInSectorAsArray, gPrimesInBlockAsArray,\
    angularDistribution, SectorRace

# cpdef gprimes(uint64_t x):
#     """Generate a GintList object of primes up to norm x."""
#     # Cython casts vector to an array.
#     return GintList(gPrimesToNorm(x), x)
#
# cpdef gprimes_sector(uint64_t x, double alpha, double beta):
#     """Generate a GintList object of primes up to norm x in given sector."""
#     return GintList(gPrimesInSector(x, alpha, beta), x, alpha, beta)
#
# cpdef gprimes_block(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy):
#     """Generate a GintList object of primes up to norm x in given sector."""
#     return GintList(gPrimesInBlock(x, y, dx, dy), x, y, dx, dy)

cpdef count_gprimes(uint64_t x):
    """Count primes up to norm x."""
    return gPrimesToNormCount(x)

cpdef count_gprimes_sector(uint64_t x, double alpha, double beta):
    """Count primes up to norm x."""
    return gPrimesInSectorCount(x, alpha, beta)

cpdef count_gprimes_block(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy):
    """Count primes up to norm x."""
    return gPrimesInBlockCount(x, y, dx, dy)

cpdef gprimes(uint64_t x):
    both = gPrimesToNormAsArray(x)
    cdef uint32_t *ptr = both.first
    cdef uint64_t size = both.second
    # Casting c++ pointer returned by gPrimesToNormArray to a memory view object
    cdef view.array primes = <cnp.uint32_t[:size]> ptr
    cdef np_primes = np.asarray(primes)
    # De-flattening array.
    np_primes = np_primes.reshape(size // 2, 2).transpose()
    return Gints(np_primes, x)

cpdef gprimes_sector(uint64_t x, double alpha, double beta):
    both = gPrimesInSectorAsArray(x, alpha, beta)
    cdef uint32_t *ptr = both.first
    cdef uint64_t size = both.second
    # Casting c++ pointer returned by gPrimesInSectorAsArray to a memory view object
    cdef view.array primes = <cnp.uint32_t[:size]> ptr
    cdef np_primes = np.asarray(primes)
    # De-flattening array.
    np_primes = np_primes.reshape(size // 2, 2).transpose()
    return Gints(np_primes, x, alpha, beta)

cpdef gprimes_block(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy):
    both = gPrimesInBlockAsArray(x, y, dx, dy)
    cdef uint32_t *ptr = both.first
    cdef uint64_t size = both.second
    # Casting c++ pointer returned by gPrimesInBlockAsArray to a memory view object
    cdef view.array primes = <cnp.uint32_t[:size]> ptr
    cdef np_primes = np.asarray(primes)
    # De-flattening array.
    np_primes = np_primes.reshape(size // 2, 2).transpose()
    return Gints(np_primes, x, y, dx, dy)

cpdef angular_dist(uint64_t x, uint32_t n_sectors):
    """Make histogram of number of Gaussian primes up to norm x in equispaced sectors."""
    data = np.array(angularDistribution(x, n_sectors))
    min = np.percentile(data, 1)
    max = np.percentile(data, 99)
    plt.hist(data, bins=30, range=(min, max))
    plt.show()
    return data



# Read https://docs.scipy.org/doc/numpy/user/basics.subclassing.html for subclasses of np.ndarray
class Gints(np.ndarray):
    """Class to wrap sieved arrays."""

    def __new__(cls, gints_np_array, x, x1=-1, x2=-1, x3=-1):
        obj = np.asarray(gints_np_array).view(cls)
        obj.x = x
        if x3 != -1:
            obj.y = x1
            obj.dx = x2
            obj.dy = x3
            obj.sieve = 'block'
        elif x2 != -1:
            obj.alpha = x1
            obj.beta = x2
            obj.sieve = 'sector'
        else:
            obj.sieve = 'to_norm'
        return obj

    def __array_finalize__(self, obj):
        pass

    def to_complex(self):
        """Convert 2d array into a 1d array of complex numbers."""
        return self[0, :] + 1j * self[1, :]

    def plot(self, full_disk=False, save=False):
        """Plot Gaussian primes with Matplotlib."""

        reals = self[0]
        imags = self[1]
        plt.subplots(figsize=(8, 8))

        if self.sieve == 'block':
            # Plotting block.
            plt.plot(reals, imags, 'ro', markersize=400 / max(self.dx, self.dy))
        else:
            # Drawing x and y axes.
            plt.axhline(0, color='red')
            plt.axvline(0, color='red')

            if self.sieve == 'to_norm' and full_disk:
                # TODO: possibly don't need to cast since changing gint struct
                reals.dtype = np.int32  # casting from unsigned to signed
                imags.dtype = np.int32

                plt.plot(reals, imags, 'bo', markersize=100 / (self.x ** .5))
                plt.plot(imags, -reals, 'bo', markersize=100 / (self.x ** .5))
                plt.plot(-reals, -imags, 'bo', markersize=100 / (self.x ** .5))
                plt.plot(-imags, reals, 'bo', markersize=100 / (self.x ** .5))

            else:
                plt.plot(reals, imags, 'bo', markersize=200 / (self.x ** .5))
                if self.sieve == 'sector':
                    for angle in [self.alpha, self.beta]:
                        plt.plot([0, np.cos(angle) * self.x ** 0.5], [0, np.sin(angle) * self.x ** 0.5], 'g-')

        if save:
            plt.savefig('sieve_visual.png')

        plt.show()



class SectorRaceWrapper:
    """Wrapper class to hold data on Gaussian prime races."""

    def __init__(self, x, alpha, beta, gamma, delta, n_bins=1000):
        if (beta - alpha) - (delta - gamma) > 0.0000001:
            raise ValueError('Unfair race; try again with two equal-sized intervals.')
        for angle in [alpha, beta, gamma, delta]:
            if angle < 0 or angle > np.pi / 4:
                raise ValueError('Stay inside the first octant.')

        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.delta = delta
        self.x = x

        self.n_bins = n_bins
        self.bin_width = self.x // self.n_bins
        self.norms = np.linspace(self.bin_width, self.x, self.n_bins)
        normalizer = lambda x: np.sqrt(x) / np.log(x)
        self.normalize = normalizer(self.norms)

        # Taking what is needed from cpp class
        race = new SectorRace(x, n_bins, alpha, beta, gamma, delta)

        s = race.getFirstSector()
        cdef uint32_t *ptr = s.first
        cdef uint64_t size = s.second
        cdef view.array primes = <cnp.uint32_t[:size]> ptr  # casting to a memory view object
        # De-flattening and chaining methods to get 2d numpy array.
        self.sector1 = np.asarray(primes).reshape(size // 2, 2).transpose()

        s = race.getSecondSector()
        ptr = s.first
        size = s.second
        primes = <cnp.uint32_t[:size]> ptr  # casting to a memory view object
        # De-flattening and chaining methods to get 2d numpy array.
        self.sector2 = np.asarray(primes).reshape(size // 2, 2).transpose()

        cdef int32_t * race_ptr = race.getNormData()
        cdef view.array norm_data = <cnp.int32_t[:n_bins]> race_ptr  # casting to a memory view object
        self.norm_data = np.asarray(norm_data)

        # Deleting cpp instance
        del race

    def plot_race(self, normalize=True):
        """Plot norms against the difference pi(sector1) - pi(sector2)."""
        plt.subplots(figsize=(8, 8))
        plt.plot(self.norms, self.norm_data, 'b-')
        if normalize:
            plt.plot(self.norms, self.normalize, 'r-')
            plt.plot(self.norms, -self.normalize, 'r-')

        plt.title('Gaussian prime race in sectors')
        plt.xlabel('$x$')
        plt.ylabel('$\pi(x; {}, {}) - \pi(x; {}, {})$'.format(self.alpha, self.beta, self.gamma, self.delta))
        plt.axhline(0, color='red')
        plt.show()

    def plot_shanks(self, bins='auto'):
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

    def plot_sectors(self, save=False):
        """Plot sectors in complex plane."""

        plt.subplots(figsize=(8, 8))
        plt.plot(self.sector1[0], self.sector1[1], 'ro', markersize=200 / (self.x ** .5))
        plt.plot(self.sector2[0], self.sector2[1], 'bo', markersize=200 / (self.x ** .5))

        for angle in [self.alpha, self.beta, self.gamma, self.delta]:
            plt.plot([0, np.cos(angle) * self.x ** 0.5], [0, np.sin(angle) * self.x ** 0.5], 'g-')

        if save:
            plt.savefig('sieve_visual.png')

        plt.show()
