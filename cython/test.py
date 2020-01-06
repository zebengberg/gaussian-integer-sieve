import gintsieve

P = gintsieve.gprimes(10000)

gintsieve.visualize(P)

print(gintsieve.gprimes_count(10000))

print(gintsieve.gprimes_region(1000, 1000, 10))

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

