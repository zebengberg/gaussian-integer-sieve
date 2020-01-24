import gintsieve
import numpy as np
import matplotlib.pyplot as plt

#gintsieve.angular_dist(10 ** 8, 10000)

gintsieve.sector_race(10 ** 9, .1, .2, .2, .3, 1000)

# p = gintsieve.gprimes_as_np(10 ** 9)
# norms = p[0] ** 2 + p[1] ** 2
# angles = np.arctan2(p[1], p[0])
# first = np.where((angles >= .1) & (angles < .2), 1, 0).cumsum()
# second = np.where((angles >= .2) & (angles < .3), 1, 0).cumsum()
# race = first - second
#
# plt.subplots(figsize=(8, 8))
# plt.plot(norms, race, 'b-')
# plt.show()





# "plt.savefig('race2', dpi=1200)"