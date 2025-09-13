import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x, y)
Z = -1 / np.sqrt(X**2 + Y**2)

# Load the saved colormap data
loaded_cmap_data = np.load("my_nipy_spectral.npy")

# Create a new colormap from loaded data
my_nipy_spectral = mcolors.LinearSegmentedColormap.from_list("my_nipy_spectral", loaded_cmap_data)

# Register it in Matplotlib for direct use
plt.register_cmap(cmap=my_nipy_spectral)

# Now you can use it like a built-in colormap!
plt.figure()
plt.contourf(X, Y, Z, levels=500, cmap="my_nipy_spectral")
plt.axis('equal')
plt.axis('square')
plt.colorbar()
plt.show()
