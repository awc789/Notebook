from main import draw_blackwell_macqueen_urn
import numpy as np
import matplotlib.pyplot as plt

# Set the default font size:
title_size = 24
axis_size = 22
legend_size = 20

if __name__ == '__main__':
    sample_times = 10
    n_draws = 10000

    for alpha in [1, 10, 100]:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for sample in range(sample_times):
            initial_counts = list()
            counts, urns = draw_blackwell_macqueen_urn(initial_counts=initial_counts, alpha=alpha, n_draws=n_draws)
            if sample != sample_times - 1:
                ax.plot(range(n_draws), urns)

        ax.plot(range(n_draws), urns)
        ax.set_title('Blackwell-MacQueen urns / alpha=%s' % alpha, size=title_size)

        ax.set_ylabel('Number of urns', size=axis_size)
        ax.set_xlabel('Number of balls', size=axis_size)
        ax.tick_params(labelsize=axis_size - 10)
        plt.legend(prop={'size': legend_size})
        plt.tight_layout()
        plt.show()