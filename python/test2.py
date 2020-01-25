import gintsieve
import numpy as np
import matplotlib.pyplot as plt

#gintsieve.angular_dist(10 ** 8, 10000)

#gintsieve.sector_race(10 ** 9, .1, .2, .2, .3, 1000)

r = gintsieve.SectorRace(1000000000, .2, .3, .4, .5, 10000000)

print(r.density())
#r.hist()
r.plot_race()