# cython: language_level=3

import matplotlib.pyplot as plt
import numpy as np
from gintsieve_externs cimport gPrimes, gPrimesCount, gPrimesSegment, gPrimesSegmentCount



# cython casts vector to a list
# see https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html#standard-library


cpdef gprimes(long x, long y=0, long z=0):
    """Generate a list of primes up to norm x or in given rectangle."""
    if z:
       return GintList(gPrimesSegment(x, y, z), x, y, z)
    else:
        return GintList(gPrimes(x), x)


cpdef gprimes_count(long x, long y=0, long z=0):
    """Count primes up to norm x or in given rectangle."""
    if z:
        return gPrimesSegmentCount(x, y, z)
    else:
        return gPrimesCount(x)


class GintList(list):
    """Class to wrap sieved lists."""
    def __init__(self, gints, x, y=0, z=0):
        list.__init__(self, gints)
        self.x = x
        if z:
            self.y = y
            self.z = z
            self.sieve = 'segment'
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

