import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from IPython.display import clear_output

# Set the default figure size:
from IPython.core.display import display, HTML

display(HTML("<style>.container { width:80% !important; }</style>"))

# Set the default font size:
title_size = 24
axis_size = 22
legend_size = 20

'----------------------------------------------------------------------------------------------------------------------'
'---------------------------------------------- GEM Distribution ------------------------------------------------------'
'----------------------------------------------------------------------------------------------------------------------'


# Using scipy.stats.beta.rvs() to generate V_i samples
# This is the process of stick-breaking
def Generate_GEM_Samples(alpha, k_trunc=10):
    """
    Generate GEM samples
    :param alpha: a positive real number
    :param k_trunc: the number of pi_i components
    :return: return a 1D array with shape (k_trunc,)
    """
    beta_dist = stats.beta(a=1, b=alpha)
    total_prob = 1.
    remain_prob = total_prob
    gem_sample = list()
    for i in range(k_trunc):
        V_i = beta_dist.rvs(1)[0]
        pi_i = V_i * remain_prob
        gem_sample.append(pi_i)
        remain_prob *= (1 - V_i)
    return np.array(gem_sample, dtype=np.float64)[:, 0]  # return a 1D array with shape (k_trunc,)


# Using the function above to obtain a Dirichlet distribution
def Dirichlet_Distribution(alpha, k_trunc, num_samples=1):
    """
    Generate a GEM distribution with k_trunc components for num_samples times
    :param alpha: a positive real number
    :param k_trunc: the number of pi_i components
    :param num_samples: the times of sampling
    :return: return a 1D array with shape (k_trunc * num_samples,)
    """
    gem_samples = list()
    for i in range(num_samples):
        sample = Generate_GEM_Samples(alpha=alpha, k_trunc=k_trunc)
        gem_samples.append(sample)
        print('Processed samples: %s/%s' % (i + 1, num_samples))
        clear_output(wait=True)
    return np.array(gem_samples, dtype=np.float64)  # return a 1D array with shape (k_trunc * num_samples,)


'----------------------------------------------------------------------------------------------------------------------'
'------------------------------------------------- Polya Urn ----------------------------------------------------------'
'----------------------------------------------------------------------------------------------------------------------'


def draw_balck_white_ball(state):
    """
    Draw a ball from the urn with the probability of each ball being drawn is proportional to its number
    :param state: a dictionary of integers
    :return: return a tuple of (draw, state)
    """
    unif = stats.uniform()
    r = unif.rvs(1)[0]
    if r < state['black'] / sum(state.values()):
        state['black'] += 1  # update the state
    else:
        state['white'] += 1  # update the state
    draw = state['black'] / sum(state.values())  # the probability of drawing a black ball
    return draw, state


def polya_urn(initial_b, initial_w, n_draws):
    """
    Generate a Polya Urn distribution with n_draws draws
    :param initial_b: the number of black balls in the urn
    :param initial_w: the number of white balls in the urn
    :param n_draws: the number of draws
    :return: return a 1D array with shape (n_draws,)
    """
    state = {'black': initial_b, 'white': initial_w}
    draws = list()
    for i in range(n_draws):
        draw, state = draw_balck_white_ball(state)
        draws.append(draw)
        # print('Processed draws: %s/%s' % (i + 1, n_draws))
        # clear_output(wait=True)
    return np.array(draws, dtype=np.float64)  # return a 1D array with shape (n_draws,)


def draw_polyas_urn(initial_b, initial_w, n_draws, num_samples):
    """
    Generate a Polya Urn distribution with n_draws draws for num_samples times
    :param initial_b: the number of black balls in the urn
    :param initial_w: the number of white balls in the urn
    :param n_draws: the number of drawsa
    :param num_samples: the times of sampling
    :return: return a 2D array with shape (n_draws, num_samples), and a 1D array with shape (num_samples,)
    """
    draws = list()
    draws_result = list()
    for i in range(num_samples):
        draw = polya_urn(initial_b, initial_w, n_draws)
        draws.append(draw)
        draws_result.append(draw[-1])
        print('Processed samples: %s/%s' % (i + 1, num_samples))
        clear_output(wait=True)
    return np.array(draws, dtype=np.float64), np.array(draws_result)
    # return a 2D array with shape (num_samples, n_draws), and a 1D array with shape (num_samples,)


'----------------------------------------------------------------------------------------------------------------------'
'------------------------------------------ Blackwell MacQueen Urn ----------------------------------------------------'
'----------------------------------------------------------------------------------------------------------------------'


def compute_urn_prob(counts, alpha):
    """
    Compute the probability of each urn
    :param counts: a list of integers, represent the number of balls in each urn
    :param alpha: a positive real number
    :return: return a 1D array with shape (len(counts),)
    """
    probs = np.array(counts)
    total_num = probs.sum() + alpha
    probs = probs / total_num
    probs = list(probs)
    probs.append(alpha / total_num)

    # compute the cumulative probability
    cum = np.zeros(len(probs))
    for i in range(len(probs)):
        cum[i] = np.sum(probs[:i + 1])
    return cum  # return a 1D array with shape (len(counts) + 1,)


def draw_new_ball_from_urn(counts, alpha):
    """
    Draw a ball from the urn with the probability of each ball being drawn is proportional to its number
    :param counts: a list of integers, represent the number of balls in each urn
    :param alpha: a positive real number
    :return: return a list of integers with length len(counts) + 1
    """
    unif = stats.uniform()
    u = unif.rvs(1)[0]
    cum = compute_urn_prob(counts, alpha)

    for i, prob_c in enumerate(cum):
        if u < prob_c:
            if i == len(cum) - 1:
                counts.append(1)
            else:
                counts[i] += 1
            break
    return counts  # return a list of integers，represent the number of balls in each urn, with length len(counts) + 1


def draw_blackwell_macqueen_urn(initial_counts, alpha, n_draws):
    """
    Generate a Blackwell MacQueen Urn distribution with n_draws draws
    :param initial_counts: a list of integers, represent the number of balls in each urn
    :param alpha: a positive real number
    :param n_draws: the number of draws
    :return: return a 1D array with shape (urns[-1],), and a 1D array with shape (n_draws,)
    """
    counts = initial_counts
    urns = list()
    for i in range(n_draws):
        counts = draw_new_ball_from_urn(counts, alpha)
        urn = len(counts)
        urns.append(urn)
        print('alpha = %s  /  Processed draws: %s/%s' % (alpha, i + 1, n_draws))
        clear_output(wait=True)
    return np.array(counts, dtype=np.float64), np.array(urns)
    # return a 1D array with shape (urns[-1],), and a 1D array with shape (n_draws,)


'----------------------------------------------------------------------------------------------------------------------'
'---------------------------------------- Chinese Restaurant Process --------------------------------------------------'
'----------------------------------------------------------------------------------------------------------------------'


def draw_new_customer_from_restaurant(initial_counts, alpha, n_customers):
    """
    Generate a Chinese Restaurant Process with n_customers
    :param initial_counts: a list of integers, represent the number of customers in each table
    :param alpha: a positive real number
    :param n_customers: the number of total customers
    :return: return two 1D array with shape (talbles[-1],), and (n_customers,)
    """
    counts = initial_counts
    talbles = list()
    for i in range(n_customers):
        # compute the cumulative probability
        prob = np.array(counts)
        total_num = prob.sum() + alpha
        prob = prob / total_num
        prob = list(prob)
        prob.append(alpha / total_num)

        cum = np.zeros(len(prob))
        for j in range(len(prob)):
            cum[j] = np.sum(prob[:j + 1])

        # draw a new customer
        unif = stats.uniform()
        u = unif.rvs(1)[0]
        for k, prob_c in enumerate(cum):
            if u < prob_c:
                if k == len(cum) - 1:
                    counts.append(1)
                else:
                    counts[k] += 1
                break

        # compute the number of tables
        table = len(counts)
        talbles.append(table)
        print('alpha = %s  /  Processed customers: %s/%s' % (alpha, i + 1, n_customers))
        clear_output(wait=True)
    return np.array(counts, dtype=np.float64), np.array(talbles)
    # return a 1D array with shape (talbles[-1],), and a 1D array with shape (n_customers,)


'----------------------------------------------------------------------------------------------------------------------'
'------------------------------------------ Dirichlet Process ---------------------------------------------------------'
'----------------------------------------------------------------------------------------------------------------------'


def dp_draw(alpha, base_measure, k_trunc):
    """
    Generate a Dirichlet Process with parameter of alpha and base_measure for k_trunc
    :param alpha: a positive real number
    :param base_measure: a function for generating a random variable which follows this base measure
    :param k_trunc: a positive integer
    :return: return a 2D array with shape (k_trunc, 2)
    """
    beta_dist = stats.beta(a=1, b=alpha)
    total_prob = 1.
    remain_prob = total_prob
    dp_sample = list()
    for i in range(k_trunc):
        V_i = beta_dist.rvs(1)[0]
        pi_i = V_i * remain_prob  # the weight of the i-th theta
        theta_i = base_measure.rvs(1)[0]  # the i-th theta
        dp_sample.append([pi_i, theta_i])
        remain_prob *= (1 - V_i)
    return np.array(dp_sample, dtype=np.float64)  # return a 2D array with shape (k_trunc, 2)


def dirichlet_process(alpha, base_measure, k_trunc, num_samples=1):
    """
    Generate a Dirichlet Process with parameter of alpha and base_measure for k_trunc
    :param alpha: a positive real number
    :param base_measure: a function for generating a random variable which follows this base measure
    :param k_trunc: a positive integer
    :param num_samples: the number of times to draw from the Dirichlet Process, default is 1
    :return: return a 3D array with shape (num_samples, k_trunc, 2)
    """
    dp_samples = np.zeros((num_samples, k_trunc, 2))
    for i in range(num_samples):
        dp_sample = dp_draw(alpha=alpha, base_measure=base_measure, k_trunc=k_trunc)
        dp_samples[i, :, :] = dp_sample
        print('Processed samples: %s/%s' % (i + 1, num_samples))
        clear_output(wait=True)

    print('DP Constructed:')
    print('\t Number of samples: %s' % dp_samples.shape[0])
    print('\t Alpha: %s' % alpha)
    print('\t Truncation value: %s' % dp_samples.shape[1])
    print('\t Bease measure: %s' % base_measure)
    return dp_samples  # return a 3D array with shape (num_samples, k_trunc, 2)













