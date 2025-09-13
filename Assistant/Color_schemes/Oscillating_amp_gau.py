import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Define the function
def func(x, y, t):
    return np.sin(t)**2 * np.exp(-x**2 - y**2)

# Create a grid
x = np.linspace(-3, 3, 100)
y = np.linspace(-3, 3, 100)
X, Y = np.meshgrid(x, y)

# Initialize the plot
fig, ax = plt.subplots()
cmap = plt.cm.jet  # Choose a colormap
norm = plt.Normalize(vmin=0, vmax=1)  # Normalize to initial function range
contour = ax.contourf(X, Y, func(X, Y, 0), levels=200, cmap=cmap, norm=norm)

# Add a colorbar
cbar = fig.colorbar(contour)
cbar.set_label("Amplitude")

# Update function for animation
def update(t):
    global contour  # Ensure we modify the global contour object
    for c in contour.collections:
        c.remove()  # Remove old contour safely
    Z = func(X, Y, t)
    contour = ax.contourf(X, Y, Z, levels=200, cmap=cmap, norm=norm)
    return contour.collections

# Create animation
ani = FuncAnimation(fig, update, frames=np.linspace(0, 2 * np.pi, 100), interval=50)

plt.show()
