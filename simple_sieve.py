###############################################################################
"""Simple sieve of Eratosthenes in the Gaussian integers.

This module implements a simple sieve of Eratosthenes in the Gaussian integers.
Methods here do not include any arithmetic statistics and only basic
visualization.  Code is written in native python as opposed to numpy or pandas.
Avoiding use of classes, and taking a purely functional approach.

Purpose of this simple sieve is two-fold.  This sieve will be used to generate
an initial sifting set of primes which will be used in more complex sieves.
Sieve will also be used as a benchmark for comparing more involved
implementations such as a donut sieve (2-dimensional wheel sieve) or wheel
sieve and a segmented sieve.  Will also compare this implementation to sieves
written in C or C++.

"""

from matplotlib import pyplot as plt


###############################################################################
def isqrt(n):
    """Integer square root based on Babylonian method."""
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x


###############################################################################
def initialize_sieve(x):
    """Create empty array corresponding to Gaussian integers with norm up to x.

    Returns G, a list of lists with boolean values initialized to True.  The
    entry G[a][b] corresponds to a Gaussian integer a + bi in the first
    quadrant.  Both indices a and b start at 0 and run up to a^2 + b^2 <= x.
    Set the elements of G corresponding to 0, 1, and i to False, because these
    are neither prime nor composite.

    """
    G = []
    for a in range(isqrt(x) + 1):
        G.append([True for b in range(isqrt(x - a**2) + 1)])
    G[0][0] = G[0][1] = G[1][0] = False
    return G


###############################################################################
def sorted_by_norm(x):
    """Gaussian integers sorted by norm with norm up to x.

    Returns a list N of pairs (a, b) corresponding to Gaussian integers a + bi
    sorted according to norm.  Here a > 0, b >= 0, and a^2 + b^2 <= x.  The
    choice of starting value of a and b makes the Gaussian integer a + bi
    unique up to action of the unit group.

    """
    N = []
    for a in range(1, isqrt(x) + 1):
        column = [(a, b) for b in range(0, isqrt(x - a ** 2) + 1)]
        N += column
    N.sort(key=lambda x: x[0] ** 2 + x[1] ** 2)
    return N


###############################################################################
def get_next_prime_index(G, N, prime_index):
    """Get the next prime index from list of sorted Gaussian integers.

    Here G is the array of Gaussian integers to be sieved, N is a list of
    sorted Gaussian integers with small norm, and prime_index is the index of
    the current prime used in sieving.

    In the rational integers, once sieving by a particular prime is over, the
    next prime is found by moving to the right along the number line.  In the
    Gaussian integers, once sieving is over, one must look outward from the
    origin to find the next prime lattice point of smallest norm.

    The list N tracks Gaussian primes of small norm, and this function finds
    the next element of N after sieving by the previous prime is complete.

    """
    is_prime = False
    while not is_prime:
        prime_index += 1
        try:
            a, b = N[prime_index]
        except IndexError:  # Beyond any pair from N.
            return False
        # Checking to see if a + bi has already been crossed from G.
        is_prime = G[a][b]
    return prime_index


###############################################################################
def cross_off_multiples(G, a, b):
    """Cross off multiples of a + bi in the sieve array G.

    Here, a > 0 and b >= 0 are the components of a Gaussian prime a + bi.  This
    prime a + bi sits above a rational prime p which can either split, ramify,
    or remain inert in Z[i].

    If b > 0, then p has degree 1, so either p = 2 or p = 1 mod 4.  In this
    case, a + bi generates an additive subgroup of order p in the additive
    group Z[i] / pZ[i].  To cross off multiples of a + bi, we will construct
    elements of this subgroup and translate them by multiples of p.

    If b = 0, then p has degree 2, so p = 3 mod 4.  In this case, a + bi = p is
    a rational integer.  There is no subgroup to consider, and we translate p
    by multiples of p.

    Remark 1: The theorem that makes this algorithm succeed when a + bi has
    degree 1 can be stated as follows.  If alpha is a Gaussian integer
    divisible by a + bi, then there exists a Gaussian integer delta and
    rational integer d such that alpha = p * delta + (a + bi) * d.

    To sketch a proof, consider the additive group (a + bi)Z[i] / pZ[i], and
    view alpha as an element of this group.  Because this group has order p, it
    is cyclic, and a + bi generates this group additively.  Therefore we have
    alpha = (a + bi) * d mod pZ[i] for some rational integer d.

    Remark 2: In the usual sieve of Eratosthenes for rational integers, because
    any multiple of p less than p^2 is divisible by some other smaller prime
    factor, we can cross off multiples of p starting at p^2.  One might ask if
    there is some version of this phenomenon in the Gaussian integers.

    If s + ti is in the subgroup generated by a split Gaussian prime a + bi,
    then N(s + ti) is a multiple of p.  If s and t are reduced mod p as in the
    loop below, then N(s + ti) < 2*p^2.  Instead, s and t could be reduced to
    lie between -p/2 and p/2.  In this case, N(s + ti) < 0.5*p^2.

    Implementing this would allow us to avoid running the loop to cross off
    the entry (s + u*p) + (t + v*p)i when both u = v = 0.  This slight
    improvement would give a negligable performance increase, and render the
    program more difficult to understand.  We avoid it.

    """

    if b:
        p = a ** 2 + b ** 2

        # Building additive subgroup generated by a + bi in Z[i]/pZ[i].
        s, t = 0, 0
        for _ in range(p):
            s += a
            s %= p
            t += b
            t %= p

            # Now translating elements of this subgroup by multiples of p.
            u = s
            for _ in range(0, len(G) - s, p):
                v = t
                for _ in range(0, len(G[u]) - t, p):
                    G[u][v] = False
                    v += p
                u += p

    else:  # So b = 0.
        u = a
        for _ in range(0, len(G) - a, a):
            v = 0
            for _ in range(0, len(G[u]), a):
                G[u][v] = False
                v += a
            u += a

    # Because we crossed off the prime a + bi, and we will want it later.
    G[a][b] = True


###############################################################################
def get_primes_after_sieve(G, sort=False):
    """Get primes after sieving."""

    P = []
    for a in range(1, len(G)):
        for b in range(len(G[a])):
            if G[a][b]:
                P.append((a, b))
    if sort:
        P.sort(key=lambda x: x[0] ** 2 + x[1] ** 2)
    return P


###############################################################################
def sieve(x, display=False, sort=False):
    """Sieve to get Gaussian primes with norm up to x."""

    if display:
        print('Initializing.')

    G = initialize_sieve(x)
    N = sorted_by_norm(isqrt(x))
    prime_index = 1  # N[0] is the pair (1, 0), which is not prime.

    if display:
        print('Starting sieve.  Each dot represents a single prime.')
        display_count = 0

    while prime_index:
        if display:
            display_count += 1
            if display_count < 100:
                print('.', end='')
            elif display_count == 100:
                print('\nNow each dot represents 10 primes.')
            elif 100 < display_count < 1000 and display_count % 10 == 0:
                print('.', end='')
            elif display_count == 1000:
                print('\nNow each dot represents 100 primes.')
            elif display_count > 1000 and display_count % 100 == 0:
                print('.', end='')

        a, b = N[prime_index]
        cross_off_multiples(G, a, b)
        prime_index = get_next_prime_index(G, N, prime_index)

    if display:
        print('\nFinishing sieve.')
        if sort:
            print('Sorting primes.')

    P = get_primes_after_sieve(G, sort)
    return P


###############################################################################
def visualize(P, x, full=False):
    """Plot Gaussian primes with matplotlib."""

    X = [p[0] for p in P]
    Y = [p[1] for p in P]
    fig, ax = plt.subplots(figsize=(15, 15))
    if full:
        neg_X = [-t for t in X]
        neg_Y = [-t for t in Y]
        plt.plot(X, Y, 'bo', markersize=200/isqrt(x))
        plt.plot(neg_X, Y, 'bo', markersize=200/isqrt(x))
        plt.plot(X, neg_Y, 'bo', markersize=200/isqrt(x))
        plt.plot(neg_X, neg_Y, 'bo', markersize=200/isqrt(x))
    else:
        plt.plot(X, Y, 'bo', markersize=300/isqrt(x))
    plt.axhline(0, color='red')
    plt.axvline(0, color='red')
    ax.set_aspect('equal')
    plt.show()




