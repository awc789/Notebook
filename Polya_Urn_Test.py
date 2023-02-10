from main import draw_polyas_urn
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# Set the default font size:
title_size = 24
axis_size = 22
legend_size = 20

if __name__ == '__main__':
    n_beta = 1000
    n_draws_polya_opt = [10, 100, 500, 1000]

    initial_b = 4
    initial_w = 1

    polya_dict = dict()
    for n_draws_polya in n_draws_polya_opt:
        _, draws = draw_polyas_urn(initial_b=initial_b, initial_w=initial_w, n_draws=n_draws_polya, num_samples=n_beta)
        polya_dict[n_draws_polya] = draws

    beta = stats.beta(a=initial_b, b=initial_w)
    range_beta = np.linspace(0, 1, 1000)
    beta_pdf = beta.pdf(range_beta)

    for n_draws_polya in polya_dict:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        draws = polya_dict[n_draws_polya]
        ax.hist(draws, 100, density=True, label='Samples of Polya Urns')
        ax.plot(range_beta, beta_pdf, label=r'$Beta(b,w)$')
        ax.set_title('Samples of independent Polya Urn \n(%s ball draws per polya)' % n_draws_polya, size=title_size)
        ax.set_ylabel(r'Probability of $\rho_{b}$', size=axis_size)
        ax.set_xlabel(r'$\rho_{b}$', size=axis_size)
        ax.tick_params(labelsize=axis_size - 10)
        plt.legend(prop={'size': legend_size})
        plt.tight_layout()