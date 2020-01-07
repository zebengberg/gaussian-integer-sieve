# cython: language_level=3

import matplotlib.pyplot as plt
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

    def visualize(self, full_disk=False):
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
            plt.plot(reals, imags, 'bo', markersize=20 / self.z)

        plt.show()