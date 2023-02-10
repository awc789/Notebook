from main import dirichlet_process
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# Set the default font size:
title_size = 24
axis_size = 22
legend_size = 20

if __name__ == '__main__':
    num_samples = 1
    base_function = stats.norm(loc=0., scale=3)
    dp_dict = dict()

    k_trunc = 10000

    for alpha in [1, 10, 100, 1000]:
        dp_samples = np.zeros((num_samples, k_trunc, 2))
        dp_samples = dirichlet_process(alpha=alpha, base_measure=base_function, k_trunc=k_trunc,
                                       num_samples=num_samples)
        dp_dict[alpha] = dp_samples

    # Plot the histogram of DP samples
    for alpha in dp_dict.keys():
        dp_samples = dp_dict[alpha]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        markeline, stemlines, baseline = ax.stem(dp_samples[0, :, 1], dp_samples[0, :, 0], '--', label='DP Draw',
                                                 use_line_collection=True)
        plt.setp(stemlines, 'color', 'b', 'linewidth', 2)
        plt.setp(baseline, 'color', 'b', 'linewidth', 0.5)

        ax.set_title(
            r'G, one draw from a DP($\alpha, H$)' + '\n' + r' $\alpha$=%s and $H=\mathcal{N}(\mu=0.,\sigma=3.)$' % alpha,
            size=title_size)
        ax.set_xlabel(r'$\Theta$ Space, support of $H$', size=axis_size)
        ax.set_ylabel(r'GEM($\alpha$=%s) probability ' % alpha + '\n' + r' assignations for each $\theta_{k}$',
                      size=axis_size)
        ax.tick_params(labelsize=axis_size - 10)
        plt.tight_layout()
        plt.legend(prop={'size': legend_size})
        plt.show()