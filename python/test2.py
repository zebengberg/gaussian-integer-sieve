import gintsieve
import numpy as np
import matplotlib.pyplot as plt

r = gintsieve.SectorRaceWrapper(10000, .1, .2, .3, .4)
print(type(r.sector1))
print(r.sector2)