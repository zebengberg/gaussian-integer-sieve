import gintsieve
import numpy as np
import matplotlib.pyplot as plt


#gintsieve.moat_explore(10000, 3.001)

x = 1200000
components = gintsieve.moat_components(x, 3.5)
first_component_size = components[0].shape[1]

fig, ax = plt.subplots(figsize=(11, 8))
for component in components:
    #color = np.random.rand(3, )  # random color
    size = component.shape[1] / first_component_size
    color = [size, 1 - size, 0]
    ax.plot(component[0], component[1], color=color, linestyle='', marker='o', markersize=400 / (x ** .5))

plt.show()