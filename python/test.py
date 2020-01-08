import gintsieve

# Print out Gaussian primes in a 10 x 10 box.
print(gintsieve.gprimes(10000, 20000, 10))

# Values taken from http://oeis.org/A091100
assert gintsieve.gprimes_count(10 ** 1) == 16
assert gintsieve.gprimes_count(10 ** 2) == 100
assert gintsieve.gprimes_count(10 ** 3) == 668
assert gintsieve.gprimes_count(10 ** 4) == 4928
assert gintsieve.gprimes_count(10 ** 5) == 38404
assert gintsieve.gprimes_count(10 ** 6) == 313752
assert gintsieve.gprimes_count(10 ** 7) == 2658344
assert gintsieve.gprimes_count(10 ** 8) == 23046512
assert gintsieve.gprimes_count(10 ** 9) == 203394764

# Graph Gaussian primes with norm up to one million.
P = gintsieve.gprimes(10 ** 6)
P.visualize(full_disk=True)
