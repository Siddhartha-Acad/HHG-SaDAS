import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

def truncate_cmap(cmap_in='jet', minval=0.0, maxval=1.0, n=100):
    cmap_in = plt.get_cmap(cmap_in)

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap_in.name, a=minval, b=maxval),
        cmap_in(np.linspace(minval, maxval, n)))

    return new_cmap
