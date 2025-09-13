import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

x0 = 0.5; a = 0.05

def sigmoid(x, x0, a):
    return 1 / (1 + np.exp(-(x - x0) / a))


def sigmoid_cmap(color_map='Blues', central_position=x0, transition_width=a, N=256):
    original_cmap = getattr(plt.cm, color_map)
    sigmoid_values = sigmoid(np.linspace(0, 1, N), central_position, transition_width)
    sigmoid_colors = original_cmap(sigmoid_values)
    return LinearSegmentedColormap.from_list("SigmoidBlues", sigmoid_colors, N=N)



if __name__ == '__main__':
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    data = np.random.rand(10, 10)
    plot1 = ax1.imshow(data, cmap='Greens')
    fig.colorbar(plot1, ax=ax1, fraction=0.046, pad=0.04)
    ax1.set_title('Normal Cmap', fontsize=15)

    plot2 = ax2.imshow(data, cmap=sigmoid_cmap(color_map='Greens', central_position=0.6, transition_width=0.01))
    fig.colorbar(plot2, ax=ax2, fraction=0.046, pad=0.04)
    ax2.set_title('Sigmoid Cmap', fontsize=15)

    plt.show()