import numpy as np
import matplotlib.pyplot as plt
# Strain = 400(-x^2 + 0.1x)
strains = np.linspace(0,0.1,29)
force_density = 400*(-strains**2 + 0.1* strains)
matrix = np.array(list(zip(strains, force_density)))
np.savetxt("force_response.csv", matrix)

# plt.plot(strains, force_density, marker = ".")
# plt.grid(True)
# plt.show()