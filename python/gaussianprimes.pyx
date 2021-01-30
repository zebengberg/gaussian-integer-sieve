# cython: language_level=3

from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t, uint64_t, int32_t
import math
import cython
import matplotlib.pyplot as plt
from cython cimport view
import numpy as np
cimport numpy as cnp
cimport gaussianprimes as gp


cdef cnp.ndarray ptr_to_np_array(pair[intptr, uint64_t] p):
  """Unpack a C++ pair holding a pointer and size into a 2D numpy array of unsigned longs."""
  cdef intptr ptr = p.first
  cdef uint64_t size = p.second
  # casting to a memory view object
  if size == 0:
    return np.array([]).reshape(2, 0)
  cdef view.array a = <cnp.int32_t[:size]> ptr
  # De-flattening and chaining methods to get 2d numpy array.
  return np.asarray(a).reshape(size // 2, 2).transpose()


# @cython.ccall is equivalent to cpdef, but allows for return annotation
# see https://github.com/cython/cython/issues/2529


cpdef count(x: int):
  """Count Gaussian primes, including associates, up to norm x.

  Args:
      x (int): Norm bound

  Returns:
      int: Count of Gaussian primes

  Raises:
      OverflowError: If x cannot be cast to uint64
  """
  return gp.gPrimesToNormCount(x)


cpdef count_sector(x: int, alpha: float, beta: float):
  """Count Gaussian primes in sector bounded by rays at alpha and beta up to norm x.

  Args:
      x (int): Norm bound
      alpha (float): Angle, in radians, of initial ray
      beta (float): Angle, in radians, of terminal ray

  Returns:
      int: Count of Gaussian primes

  Raises:
      OverflowError: If x cannot be cast to uint64
      NotImplementedError: If alpha and beta are not in the interval [0, p/2)
  """
  if min(alpha, beta) < 0 or max(alpha, beta) > math.pi / 2:
    raise NotImplementedError('Only implemented for alpha >= 0 and beta < pi/2.')
  # SectorSieve automatically assigns alpha to smaller and beta to larger
  return gp.gPrimesInSectorCount(x, alpha, beta)


cpdef count_block(x: int, y: int, dx: int, dy: int):
  """Count Gaussian primes in the block [x, x + dx) x [y, y + dy).

  Args:
      x (int): Real coordinate of lower left corner
      y (int): Imaginery coordinate of lower left corner
      dx (int): Block width
      dy (int): Block height

  Returns:
      int: Count of Gaussian primes

  Raises:
      OverflowError: If x, y, dx, or dy cannot be cast to uint32
  """
  return gp.gPrimesInBlockCount(x, y, dx, dy)


cpdef gprimes(x: int):
  """Return Gaussian primes in first quadrant up to norm x.

  Args:
      x (int): Norm bound

  Returns:
      Gints: Array of Gaussian primes

  Raises:
      OverflowError: If x cannot be cast to uint64
  """
  p = gp.gPrimesToNormAsArray(x)
  np_primes = ptr_to_np_array(p)
  return Gints(np_primes, x)


cpdef gprimes_sector(x: int, alpha: float, beta: float):
  """Return Gaussian primes in sector bounded by rays at alpha and beta up to norm x.

  Args:
      x (int): Norm bound
      alpha (float): Angle, in radians, of initial ray
      beta (float): Angle, in radians, of terminal ray

  Returns:
      Gints: Array of Gaussian primes

  Raises:
      OverflowError: If x cannot be cast to uint64
      NotImplementedError: If alpha and beta are not in the interval [0, p/2)
  """
  if min(alpha, beta) < 0 or max(alpha, beta) > math.pi / 2:
    raise NotImplementedError('Only implemented for alpha >= 0 and beta < pi/2.')
  # SectorSieve automatically assigns alpha to smaller and beta to larger
  p = gp.gPrimesInSectorAsArray(x, alpha, beta)
  np_primes = ptr_to_np_array(p)
  return Gints(np_primes, x, alpha, beta)


cpdef gprimes_block(x: int, y: int, dx: int, dy: int):
  """Return Gaussian primes in the block [x, x + dx) x [y, y + dy).

  Args:
      x (int): Real coordinate of lower left corner
      y (int): Imaginery coordinate of lower left corner
      dx (int): Block width
      dy (int): Block height

  Returns:
      Gints: Array of Gaussian primes

  Raises:
      OverflowError: If x, y, dx, or dy cannot be cast to uint32
  """
  p = gp.gPrimesInBlockAsArray(x, y, dx, dy)
  np_primes = ptr_to_np_array(p)
  return Gints(np_primes, x, y, dx, dy)


cpdef angular_dist(x: int, n: int, ignore_outliers: bool = True):
  """Create histogram of Gaussian primes up to norm x in n equal-spaced sectors.

  Args:
      x (int): Norm bound
      n (int): Number of bins
      ignore_outliers (bool): Ignore histogram outliers, default True

  Returns:
      np.array: Array of histogram counts

  Raises:
      OverflowError: If x cannot be cast to uint64 or n cannot be cast to uint32
  """
  data = np.array(gp.angularDistribution(x, n))
  m0 = np.percentile(data, 1)
  m1 = np.percentile(data, 99)
  plt.subplots(figsize=(12, 8))
  plt.title('Histogram of Gaussian primes up to norm $x$ in $n$ equal-spaced sectors', fontsize=15)
  plt.xlabel('Number of Gaussian primes per sector', fontsize=15)
  plt.ylabel('Number of sectors', fontsize=15)
  if ignore_outliers:
    plt.hist(data, bins=50, range=(m0, m1))
  else:
    plt.hist(data, bins=50)
  plt.show()
  return data


cpdef moat_component(jump_size: float):
  """Calculate the connected component of the Gaussian moat graph in the first octant.

  Args:
      jump_size (float): Two Gaussian primes are adjacent iff they have distance <= jump_size

  Returns:
      Gints: Array of Gaussian primes

  Raises:
      OverflowError: If jump_size > 5
  """
  if jump_size > 5:
    raise OverflowError('Cannot handle jump_size > 5.')

  p = gp.moatMainComponent(jump_size)
  np_primes = ptr_to_np_array(p)
  # In cython_bindings.cpp, appending the largest element to the end of the array
  # Here, getting its value so we can pass it to the Gints class
  x = np_primes[0][-1] ** 2 + np_primes[1][-1] ** 2
  # removing that final element with some slicing
  np_primes = np_primes[:, :-1]
  return Gints(np_primes, x)


cpdef moat_components_to_norm(jump_size: float, x: int):
  """Calculate all connected component of the Gaussian moat graph in the first octant up to norm x.

  Args:
      jump_size (float): Two Gaussian primes are adjacent iff they have distance <= jump_size
      x (int): Norm bound

  Returns:
      list[Gints]: List of arrays of Gaussian primes

  Raises:
      OverflowError: If x cannot be cast to uint64
      OverflowError: If jump_size > 5
  """
  if jump_size > 5:
    raise OverflowError('Cannot handle jump_size > 5.')

  # Cython compiler gets confused if this isn't explicitly typed
  cdef vector[pair[intptr, uint64_t]] vector_of_ptrs = gp.moatComponentsToNorm(jump_size, x)
  components = []
  for i in range(vector_of_ptrs.size()):
    p = vector_of_ptrs[i]
    np_primes = ptr_to_np_array(p)
    components.append(np_primes)
  return components


cpdef moat_components_in_block(jump_size: float, x: int, y: int, dx: int, dy: int, ignore_edges: bool = True):
  """Calculate all connected component of the Gaussian moat graph in the block [x, x + dx) x [y, y + dy).

  Args:
      jump_size (float): Two Gaussian primes are adjacent iff they have distance <= jump_size
      x (int): Real coordinate of lower left corner
      y (int): Imaginery coordinate of lower left corner
      dx (int): Block width
      dy (int): Block height
      ignore_edges (bool): Ignore graph edges and return only vertices

  Returns:
      list | tuple[list, list]: Vertices and edges of graph.

  Raises:
      OverflowError: If x, y, dx, dy cannot be cast to uint32
  """
  # Cython compiler gets confused if this isn't explicitly typed
  cdef vector[pair[intptr, uint64_t]] vector_of_ptrs = gp.moatComponentsInBlock(jump_size, x, y, dx, dy)
  components = []
  for i in range(vector_of_ptrs.size()):
    p = vector_of_ptrs[i]
    np_primes = ptr_to_np_array(p)
    components.append(np_primes)

  # Crudely passed edge data as last component; now crudely dealing with it
  edges = components.pop()
  if ignore_edges:
    return components

  edges = edges.transpose()
  segments = []
  for i in range(edges.shape[0] // 2):
    segment = np.stack((edges[2 * i], edges[2 * i + 1]))
    segments.append(segment)

  return components, segments


class Race:
  """Wrapper class to hold data from Gaussian prime races.

  Args:
      x (int): Sieve both sectors to find Gaussian primes with norm up to x.
      alpha (float): Initial angle for first sector.
      beta (float): Terminal angle for first sector.
      gamma (float): Initial angle for second sector.
      delta (float): Terminal angle for second sector.
      n_bins (int): Number of histogram bins. Defaults to 1000.
  """

  def __init__(self, x: int, alpha: float, beta: float, gamma: float, delta: float, n_bins=1000):
    if alpha > beta or gamma > delta:
      raise ValueError('The four angle measures must be increasing.')
    if abs(abs(beta - alpha) - abs(delta - gamma)) > 0.0000001:
      raise ValueError('Unfair race; try again with two equal-sized sectors.')
    for angle in [alpha, beta, gamma, delta]:
      if angle < 0 or angle > np.pi / 4:
        raise ValueError('Stay inside the first octant.')

    self.x = x
    self.alpha = alpha
    self.beta = beta
    self.gamma = gamma
    self.delta = delta

    self.n_bins = n_bins
    self.bin_width = self.x // self.n_bins
    self.norms = np.linspace(self.bin_width, self.x, self.n_bins)
    def normalizer(x): return np.sqrt(x) / np.log(x)
    self.normalize = normalizer(self.norms)

    # Taking what is needed from cpp class
    race = new gp.SectorRace(x, n_bins, alpha, beta, gamma, delta)

    s = race.getFirstSector()
    self.sector1 = ptr_to_np_array(s)

    s = race.getSecondSector()
    self.sector2 = ptr_to_np_array(s)

    cdef int32_t * race_ptr = race.getNormData()
    # casting to a memory view object
    cdef view.array norm_data = <cnp.int32_t[:n_bins] > race_ptr
    self.norm_data = np.asarray(norm_data)

    # Deleting cpp instance; not sure if this is useful.
    del race

  def plot_race(self, normalize=True):
    """Plot norms against the difference pi(sector1) - pi(sector2)."""
    plt.subplots(figsize=(12, 8))
    plt.plot(self.norms, self.norm_data, 'b-')
    if normalize:
      plt.plot(self.norms, self.normalize, 'r-')
      plt.plot(self.norms, -self.normalize, 'r-')

    plt.title('Gaussian prime race in sectors', fontsize=15)
    plt.xlabel('$x$', fontsize=15)
    plt.ylabel('$\pi(x; {:.2f}, {:.2f}) - \pi(x; {:.2f}, {:.2f})$'.format(self.alpha, self.beta, self.gamma, self.delta),
               fontsize=15)
    plt.axhline(0, color='red')
    plt.show()

  def plot_shanks(self, bins='auto'):
    """Plot Shanks-style histogram of race progress."""
    plt.subplots(figsize=(12, 8))
    plt.title('Shanks histogram of Gaussian prime race in sectors', fontsize=15)
    plt.xlabel(r'$\frac{\pi_1(x) - \pi_2(x)}{\sqrt{x} / \log x}$', fontsize=25)
    plt.ylabel('Proportion of occurrences', fontsize=15)
    shanks = self.norm_data / self.normalize
    plt.hist(shanks, bins=bins, density=True)
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

    plt.subplots(figsize=(12, 8))
    plt.plot(self.sector1[0], self.sector1[1],
             'ro', markersize=100 / (self.x ** .5))
    plt.plot(self.sector2[0], self.sector2[1],
             'bo', markersize=100 / (self.x ** .5))

    for angle in [self.alpha, self.beta, self.gamma, self.delta]:
      plt.plot([0, np.cos(angle) * self.x ** 0.5],
               [0, np.sin(angle) * self.x ** 0.5], 'g-')

    if save:
      plt.savefig('sieve_visual.png')

    plt.show()



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
    """Convert to a 1d np.array of complex numbers."""
    return np.asarray(self[0, :] + 1j * self[1, :])

  def norms(self):
    """Convert to 1d np.array of norms."""
    return np.asarray(self[0, :] ** 2 + self[1, :] ** 2)

  def __mod__(self, z):
    """Return the remainder upon division."""
    if not isinstance(z, complex):
      raise TypeError(
          'The modulus should be a python complex number such as 3 + 4j.')
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
      # Scaling axes
      m = max(self.dx, self.dy)
      plt.subplots(figsize=(12 * self.dx / m, 12 * self.dy / m))

      # Plotting block.
      plt.plot(reals, imags, 'ro', markersize=400 / m)

    else:
      # Scaling axes
      mx = int(np.max(self[0]))
      my = int(np.max(self[1]))
      m = max(mx, my)
      plt.subplots(figsize=(12 * mx / m, 12 * my / m))
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
            plt.plot([0, np.cos(angle) * self.x ** 0.5],
                     [0, np.sin(angle) * self.x ** 0.5], 'g-')

    if save:
      plt.savefig('sieve_visual.png')

    plt.show()