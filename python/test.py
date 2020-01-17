import gintsieve

# Print out Gaussian primes in a 10 x 20 box.
print(gintsieve.gprimes_block(10000, 20000, 10, 20))

# Print out Gaussian primes in a tiny sector
print(gintsieve.gprimes_sector(10 ** 5, .11, .12))

# Values taken from http://oeis.org/A091100
assert gintsieve.count_gprimes(10 ** 1) == 16
assert gintsieve.count_gprimes(10 ** 2) == 100
assert gintsieve.count_gprimes(10 ** 3) == 668
assert gintsieve.count_gprimes(10 ** 4) == 4928
assert gintsieve.count_gprimes(10 ** 5) == 38404
assert gintsieve.count_gprimes(10 ** 6) == 313752
assert gintsieve.count_gprimes(10 ** 7) == 2658344
assert gintsieve.count_gprimes(10 ** 8) == 23046512
assert gintsieve.count_gprimes(10 ** 9) == 203394764

# Graph Gaussian primes with norm up to one million.
P = gintsieve.gprimes(10 ** 6)
P.visualize(full_disk=True)
