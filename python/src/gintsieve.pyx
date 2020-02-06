# cython: language_level=3

import numpy as np
cimport numpy as cnp
from cython cimport view
import matplotlib.pyplot as plt
from libc.stdint cimport uint32_t, uint64_t, int32_t
from libcpp.pair cimport pair

from gintsieve_externs cimport gPrimesToNorm, gPrimesInSector, gPrimesInBlock,\
    gPrimesToNormCount, gPrimesInSectorCount, gPrimesInBlockCount,\
    gPrimesToNormAsArray, gPrimesInSectorAsArray, gPrimesInBlockAsArray,\
    angularDistribution, SectorRace, OctantMoat




ctypedef int32_t * intptr
cdef cnp.ndarray ptr_to_np_array(pair[intptr, uint64_t] p):
    """Unpack a C++ pair holding a pointer and size into a 2D numpy array of unsigned longs."""
    cdef intptr ptr = p.first
    cdef uint64_t size = p.second
    cdef view.array a = <cnp.int32_t[:size]> ptr  # casting to a memory view object
    # De-flattening and chaining methods to get 2d numpy array.
    return np.asarray(a).reshape(size // 2, 2).transpose()





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

# TODO: consider removing above functions, test everything, cnp vs np and cdef, cpdef, etc. be explicit

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
    p = gPrimesToNormAsArray(x)
    np_primes = ptr_to_np_array(p)
    return Gints(np_primes, x)

cpdef gprimes_sector(uint64_t x, double alpha, double beta):
    p = gPrimesInSectorAsArray(x, alpha, beta)
    np_primes = ptr_to_np_array(p)
    return Gints(np_primes, x, alpha, beta)

cpdef gprimes_block(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy):
    p = gPrimesInBlockAsArray(x, y, dx, dy)
    np_primes = ptr_to_np_array(p)
    return Gints(np_primes, x, y, dx, dy)

cpdef angular_dist(uint64_t x, uint32_t n_sectors):
    """Make histogram of number of Gaussian primes up to norm x in equispaced sectors."""
    data = np.array(angularDistribution(x, n_sectors))
    min = np.percentile(data, 1)
    max = np.percentile(data, 99)
    plt.hist(data, bins=30, range=(min, max))
    plt.show()
    return data



# Read https://docs.scipy.org/doc/numpy/user/basics.subclassing.html for subclassing np.ndarray
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
                # Rotating around the four quadrants.
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
        self.sector1 = ptr_to_np_array(s)

        s = race.getSecondSector()
        self.sector2 = ptr_to_np_array(s)

        cdef int32_t * race_ptr = race.getNormData()
        cdef view.array norm_data = <cnp.int32_t[:n_bins]> race_ptr  # casting to a memory view object
        self.norm_data = np.asarray(norm_data)

        # Deleting cpp instance; not sure if this is useful.
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


cpdef moat_explore(x, jump_size):
    """Wrapper function to hold Gaussian moat explorations."""

    # Taking what is needed from cpp class
    moat = new OctantMoat(x, jump_size)
    moat.explore()

    p = moat.getExplored()
    explored = ptr_to_np_array(p)

    p = moat.getUnexplored()
    unexplored = ptr_to_np_array(p)

    plt.subplots(figsize=(11, 8))  # ratio should be sqrt(2) : 1
    plt.plot(explored[0], explored[1], 'ro', markersize=200 / (x ** .5))
    plt.plot(unexplored[0], unexplored[1], 'bo', markersize=200 / (x ** .5))
    plt.show()


cpdef moat_components(x, jump_size):
    # Taking what is needed from cpp class
    moat = new OctantMoat(x, jump_size)

    # Super slow conversion from vector of vector of pairs of ints to
    # list of list of tuples of python numbers.
    components = moat.getAllComponents()
    return components

