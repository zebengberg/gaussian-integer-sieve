import gintsieve
import numpy as np
import matplotlib.pyplot as plt


r = gintsieve.SectorRaceWrapper(1000000, .1, .3, .4, .6)
r.plot_sectors()
r.plot_race()
print(r.density())
r.plot_shanks()