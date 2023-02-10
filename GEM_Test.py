from main import Dirichlet_Distribution
import numpy as np
import matplotlib.pyplot as plt

# Set the default font size:
title_size = 24
axis_size = 22
legend_size = 20

if __name__ == '__main__':
    num_samples = 1000
    GEM_dict = dict()

    k_trunc = 1000
    alpha = 100
    gem_samples_big_alpha = Dirichlet_Distribution(alpha, k_trunc, num_samples)
    GEM_dict[alpha] = gem_samples_big_alpha

    k_trunc = 100
    alpha = 3
    gem_samples_small_alpha = Dirichlet_Distribution(alpha, k_trunc, num_samples)
    GEM_dict[alpha] = gem_samples_small_alpha

    # Plot the histogram of GEM samples
    dimensions = 20
    for alpha in GEM_dict.keys():
        gem_samples = GEM_dict[alpha]
        x = list(range(dimensions + 2))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.violinplot(gem_samples[:, :dimensions + 1])
        ax.set_title(r'GEM($\alpha$=%s) dimension distributions' % alpha, size=title_size)
        ax.set_ylabel('Distribution for \n each dimension', size=axis_size)
        ax.set_xlabel(r'Dimension $i$ of GEM / Cluster weight', size=axis_size)
        ax.set_xticks(x)
        plt.tight_layout()
