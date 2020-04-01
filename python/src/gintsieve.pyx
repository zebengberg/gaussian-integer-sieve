# cython: language_level=3

import numpy as np
cimport numpy as cnp
from cython cimport view
import matplotlib.pyplot as plt
from libc.stdint cimport uint32_t, uint64_t, int32_t
from libcpp.pair cimport pair
from libcpp.vector cimport vector

from gintsieve_externs cimport gPrimesToNormCount, \
                               gPrimesInSectorCount, \
                               gPrimesInBlockCount, \
                               gPrimesToNormAsArray, \
                               gPrimesInSectorAsArray, \
                               gPrimesInBlockAsArray, \
                               angularDistribution, \
                               SectorRace, \
                               moatMainComponent, \
                               moatComponentsToNorm, \
                               moatComponentsInBlock




ctypedef int32_t * intptr
cdef cnp.ndarray ptr_to_np_array(pair[intptr, uint64_t] p):
    """Unpack a C++ pair holding a pointer and size into a 2D numpy array of unsigned longs."""
    cdef intptr ptr = p.first
    cdef uint64_t size = p.second
    cdef view.array a = <cnp.int32_t[:size]> ptr  # casting to a memory view object
    # De-flattening and chaining methods to get 2d numpy array.
    return np.asarray(a).reshape(size // 2, 2).transpose()

cpdef count_gprimes(uint64_t x):
    """Count primes in first octant up to norm x."""
    return gPrimesToNormCount(x)

cpdef count_gprimes_in_sector(uint64_t x, double alpha, double beta):
    """Count primes in sector between rays at alpha and beta up to norm x."""
    return gPrimesInSectorCount(x, alpha, beta)

cpdef count_gprimes_in_block(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy):
    """Count primes in the block [x, x + dx) x [y, y + dy)."""
    return gPrimesInBlockCount(x, y, dx, dy)

cpdef gprimes(uint64_t x):
    """Get primes in first octant up to norm x."""
    if x < 2:
        raise ValueError('No primes with norm under 2!')
    p = gPrimesToNormAsArray(x)
    np_primes = ptr_to_np_array(p)
    return Gints(np_primes, x)

cpdef gprimes_in_sector(uint64_t x, double alpha, double beta):
    """Get primes in sector between rays at alpha and beta up to norm x."""
    if x < 2:
        raise ValueError('No primes with norm under 2!')
    p = gPrimesInSectorAsArray(x, alpha, beta)
    np_primes = ptr_to_np_array(p)
    return Gints(np_primes, x, alpha, beta)

cpdef gprimes_in_block(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy):
    """Count primes in the block [x, x + dx) x [y, y + dy)."""
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
    """Class to wrap np arrays returned by sieve calls."""

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

    def to_complex(self):
        """Convert 2d np array into a 1d np array of complex numbers."""
        return self[0, :] + 1j * self[1, :]

    def __mod__(self, z):
        """Return the remainder upon division."""
        if not isinstance(z, complex):
            raise TypeError('The modulus should be a python complex number such as 3 + 4j.')
        a = round(z.real)
        b = round(z.imag)
        if a != z.real or b != z.imag:
            raise ValueError('The modulus should have integer coordinates!')
        if a == 0 and b == 0:
            raise ZeroDivisionError('The modulus should not be 0!')

        # Using division algorithm in Z[i]
        # Find K Conrad paper for nice exposition
        norm = a * a + b * b
        real_quot = np.round((a * self[0] + b * self[1]) / norm)
        imag_quot = np.round((a * self[1] - b * self[0]) / norm)
        real_mod = self[0] - a * real_quot + b * imag_quot
        imag_mod = self[1] - a * imag_quot - b * real_quot

        return np.stack((real_mod.astype(np.int32), imag_mod.astype(np.int32)))


    def plot(self, full_disk=False, save=False):
        """Plot Gaussian primes with Matplotlib."""

        reals = self[0]
        imags = self[1]


        if self.sieve == 'block':
            # Plotting block.
            m = max(self.dx, self.dy)
            plt.subplots(figsize=(12 * self.dx / m, 12 * self.dy / m))
            plt.plot(reals, imags, 'ro', markersize=400 / m)
        else:
            plt.subplots(figsize=(12, 12))
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
                plt.plot(reals, imags, 'bo', markersize=300 / (self.x ** .5))
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


# Several functions for exploring the Gaussian moat graph

cpdef moat_main_component(double jump_size):
    """Calculate the main connected component of the Gaussian moat graph in the first octant starting at origin."""
    p = moatMainComponent(jump_size)
    np_primes = ptr_to_np_array(p)
    # In cython_bindings.cpp, appending the largest element to the end of the array
    # Here, getting its value so we can pass it to the Gints class
    x = np_primes[-1][0] ** 2 + np_primes[-1][1] ** 2
    np_primes = np_primes[:-1]  # removing that final element
    return Gints(np_primes, x)


cpdef moat_components_to_norm(double jump_size, uint64_t x):
    """Calculate all connected components of the Gaussian moat graph in the first octant up to norm x."""
    # Cython compiler gets confused if this isn't explicitly typed
    cdef vector[pair[intptr, uint64_t]] vector_of_ptrs = moatComponentsToNorm(jump_size, x)
    components = []
    for i in range(vector_of_ptrs.size()):
        p = vector_of_ptrs[i]
        np_primes = ptr_to_np_array(p)
        components.append(np_primes)
    return components


cpdef moat_components_in_block(double jump_size, int32_t x, int32_t y, int32_t dx, int32_t dy):
    """Calculate all connected components of the Guassian moat graph in the block [x, x + dx) x [y, y + dy)."""
    # Cython compiler gets confused if this isn't explicitly typed
    cdef vector[pair[intptr, uint64_t]] vector_of_ptrs = moatComponentsInBlock(jump_size, x, y, dx, dy)
    components = []
    for i in range(vector_of_ptrs.size()):
        p = vector_of_ptrs[i]
        np_primes = ptr_to_np_array(p)
        components.append(np_primes)
    return components