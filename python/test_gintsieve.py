import gintsieve
print('\nTesting the module just built...')


def test_count_gprimes():
  # Values taken from http://oeis.org/A091100
  assert gintsieve.count_gprimes(10 ** 0) == 0
  assert gintsieve.count_gprimes(10 ** 1) == 16
  assert gintsieve.count_gprimes(10 ** 2) == 100
  assert gintsieve.count_gprimes(10 ** 3) == 668
  assert gintsieve.count_gprimes(10 ** 4) == 4928
  assert gintsieve.count_gprimes(10 ** 5) == 38404
  assert gintsieve.count_gprimes(10 ** 6) == 313752
  assert gintsieve.count_gprimes(10 ** 7) == 2658344
  assert gintsieve.count_gprimes(10 ** 8) == 23046512
  assert gintsieve.count_gprimes(10 ** 9) == 203394764
