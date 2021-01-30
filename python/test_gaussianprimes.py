"""Test gaussianprimes module."""

import numpy as np
from sympy import isprime
import gaussianprimes as gp


def test_count_gprimes():
  """Test gp.counts by comparing to http://oeis.org/A091100 ."""
  assert gp.count(10 ** 0) == 0
  assert gp.count(10 ** 1) == 16
  assert gp.count(10 ** 2) == 100
  assert gp.count(10 ** 3) == 668
  assert gp.count(10 ** 4) == 4928
  assert gp.count(10 ** 5) == 38404
  assert gp.count(10 ** 6) == 313752
  assert gp.count(10 ** 7) == 2658344
  assert gp.count(10 ** 8) == 23046512
  assert gp.count(10 ** 9) == 203394764

  assert gp.count(0) == 0
  assert gp.count(1) == 0
  assert gp.count(2) == 4
  assert gp.count(3) == 4
  assert gp.count(4) == 4
  assert gp.count(5) == 12
  assert gp.count(6) == 12
  assert gp.count(7) == 12
  assert gp.count(8) == 12
  assert gp.count(9) == 16

  try:
    gp.count(-1)
    raise ValueError
  except OverflowError:
    pass

  try:
    gp.count(2 ** 65)
    raise ValueError
  except OverflowError:
    pass


def test_count_block():
  """Test gp.count_block by considering page 10 on
  https://www.maa.org/sites/default/files/pdf/upload_library/22/Chauvenet/Gethner.pdf"""
  assert gp.count_block(20785207 - 10, 0, 20, 10) == 1


def verify_splitting(gints):
  """Verify that the primes obtained sit above a rational prime."""
  for p, n in zip(gints.transpose(), gints.norms()):
    if p[1]:
      assert isprime(n)
    else:
      assert isprime(p[0])


def test_gprimes():
  """Test gprimes."""
  g = gp.gprimes(10000)
  verify_splitting(g)

  g = gp.gprimes(0)
  assert g.shape == (2, 0)

  g = gp.gprimes(10)
  g = np.asarray(g)
  a = np.array([[1, 1], [2, 1], [1, 2], [3, 0]])
  a = a.transpose()
  assert (g == a).all()

  try:
    gp.gprimes(-1)
    raise ValueError
  except OverflowError:
    pass

  try:
    gp.gprimes(2 ** 65)
    raise ValueError
  except OverflowError:
    pass


def test_gprimes_block():
  """Test gprimes_block."""
  g = gp.gprimes_block(10000, 20000, 300, 400)
  verify_splitting(g)

  # from https://oeis.org/A002496/list
  p = [2, 5, 17, 37, 101, 197, 257, 401, 577, 677, 1297, 1601, 2917,
       3137, 4357, 5477, 7057, 8101, 8837, 12101, 13457, 14401,
       15377, 15877, 16901, 17957, 21317, 22501, 24337, 25601,
       28901, 30977, 32401, 33857, 41617, 42437, 44101, 50177]
  g = gp.gprimes_block(0, 1, 225, 1)
  n = g.norms().tolist()
  assert p == n


def test_gprimes_sector():
  """Test gprimes_sector."""
  g = gp.gprimes_sector(100000, 0.1, 0.2)
  verify_splitting(g)


def test_moat():
  """Test main moat function."""
  # values from https://www.maa.org/sites/default/files/pdf/upload_library/22/Chauvenet/Gethner.pdf

  m = gp.moat_component(3)
  assert m.shape == (2, 380)

  m = gp.moat_component(3.2)
  assert m.shape == (2, 31221)

  m = gp.moat_component(4)
  assert m.shape == (2, 347638)

  m = gp.moat_component(4.3)
  assert m.shape == (2, 2386129)


def test_readme_examples():
  """Test examples in readme."""
  print(gp.gprimes(50))
  print(gp.count(3141592653))
  p = gp.gprimes_block(123456, 67890, 100, 100)
  p.plot()
  p = gp.gprimes_sector(10000, 0, 3.14159 / 8)
  p.plot()
  p = gp.gprimes(1000)
  p.plot(full_disk=True)


if __name__ == '__main__':
  test_count_gprimes()
  test_count_block()
  test_gprimes()
  test_gprimes_block()
  test_gprimes_sector()
  test_moat()
  test_readme_examples()
