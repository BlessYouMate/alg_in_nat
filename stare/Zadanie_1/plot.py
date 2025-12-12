# plot_histograms.py
import numpy as np
import matplotlib.pyplot as plt
import math

def plot_for_file(datafile, N, bins=20, m=0.0, sigma=1.0, name=""):
    data = np.loadtxt(datafile)            # plik geneated by C++
    minVal = m - 3.5*sigma
    maxVal = m + 3.5*sigma
    edges = np.linspace(minVal, maxVal, bins+1)
    counts, _ = np.histogram(data, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2
    width = edges[1] - edges[0]

    plt.bar(centers, counts, width=width, edgecolor='black', alpha=0.6, label='data counts')

    x = np.linspace(minVal, maxVal, 400)
    pdf = (1.0 / (math.sqrt(2*math.pi) * sigma)) * np.exp(-0.5 * ((x - m)/sigma)**2)
    expected_counts = pdf * N * width   # E[count in bin] ≈ N * f(center) * bin_width
    plt.plot(x, expected_counts, linewidth=2, label='theoretical (expected counts)')

    plt.title(f'Histogram {name} (N={N})')
    plt.xlabel('x')
    plt.ylabel('counts')
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_uniform_for_file(datafile, N, bins=20, name=""):
    data = np.loadtxt(datafile)

    minVal = 0.0
    maxVal = 1.0
    edges = np.linspace(minVal, maxVal, bins + 1)
    counts, _ = np.histogram(data, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2
    width = edges[1] - edges[0]

    plt.bar(centers, counts, width=width, edgecolor='black', alpha=0.6, label='data counts')

    x = np.linspace(minVal, maxVal, 400)
    pdf = np.ones_like(x)  # f(x) = 1 dla x∈[0,1]
    expected_counts = pdf * N * width
    plt.plot(x, expected_counts, 'r-', linewidth=2, label='theoretical (expected counts)')

    plt.title(f'Uniform Histogram {name} (N={N})')
    plt.xlabel('x')
    plt.ylabel('counts')
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_for_file("data_10k_2_data.txt", 10000, bins=20, name="10k")
    plot_for_file("data_100k_2_data.txt", 100000, bins=20, name="100k")
    plot_for_file("data_1M_2_data.txt", 1000000, bins=20, name="1M")
    plot_uniform_for_file("data_10k_1_data.txt", 10000, bins=20, name="10k")
    plot_uniform_for_file("data_100k_1_data.txt", 100000, bins=20, name="100k")
    plot_uniform_for_file("data_1M_1_data.txt", 1000000, bins=20, name="1M")
