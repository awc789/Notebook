library(gtools)

# the value of the variables for the subjects is arranged in the matrix x
x <- matrix(c(
  0, 0, 0,
  0, 0, 0,
  0, 0, 0,
  0, 0, 0,
  0, 0, 0,
  1, 1, 1,
  1, 1, 1,
  1, 1, 1,
  1, 1, 1,
  1, 1, 1
), nrow = 10, byrow = TRUE)

iter <- 1000
burnin <- 5000
right_limit_alphadp <- 6 # uniform prior for alphadp U(0.0001,right_limit_alphadp)

lambdap <- 0.5 # hyperparameter for the prior of phi / Dirichlet distribution

N <- 10 # number of observed subjects
C <- 10 # number of clusters
P <- 3 # number of variables

# iter x N   --> total int(iter) samples of z
z <- matrix(rep(0, iter * N), nrow = iter, ncol = N) # z_i == c_i ,which is the clusters of the subjects
z_current <- matrix(rep(0, N), nrow = 1, ncol = N) # 1 x N   --> current value of z
z_new <- matrix(rep(0, N), nrow = 1, ncol = N) # 1 x N   --> new value of z

no_of_levels <- 2 # number of levels for the variables
# for the StBreaking we sample the whole phi vector
# phi
phi <- array(0, dim = c(iter, C, P, no_of_levels)) # shape: iter x C x P x noOfLevels --> iter x C for (P x noOfLevels) matrices
phi_current <- array(0, dim = c(C, P, no_of_levels)) # shape: C x P x noOfLevels --> C x P for (noOfLevels) matrices
# psi
psi <- array(0, dim = c(iter, C)) # shape: iter x C matrix
psi_current <- matrix(rep(0, C), nrow = 1, ncol = C) # shape: 1 x C matrix
# alpha == 1
alphadp <- rep(1, iter) # 1 for iter times
alphadp_current <- 1

# random group allocation to start
for (i_sub in 1:N) {
  # random generation number from 1 to C with equal probability, and assign it to z_current for N times
  #      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
  # [1,]    5    8    1   10    4   10    6    2    3     4
  z_current[i_sub] <- which(as.vector(rmultinom(1, 1, rep(1 / C, C)) == 1))
}
# the standard deviation for the proposal distribution for alphadp
st_dev_for_alphadp_prop <- 0.2

for (i_iter in 2:(iter + burnin)) { # i_iter from 2 to (iter + burnin)
  print(i_iter)
  # sample \phi
  for (i_cluster in 1:C) {
    for (ivar in 1:P) {
      number_of_successes <- sum(x[z_current == i_cluster, ivar] == 1)
      phi_current[i_cluster, ivar, ] <- rdirichlet(1, c(lambdap + number_of_successes, sum(z_current == i_cluster) + lambdap - number_of_successes))
    }
  }

  # sample \psi
  for_psi_probs <- sum(z_current == 1)
  for (i_psi_probs in 2:C) {
    for_psi_probs <- c(for_psi_probs, sum(z_current == i_psi_probs))
  }
  psi_current <- rdirichlet(1, for_psi_probs + alphadp_current)
  # if any of the psi_current is zero, then we add a small number to it
  if ((sum(psi_current < (10^(-10)))) > 0) {
    sum_of_psi_zero <- sum(psi_current < (10^(-10)))
    psi_current[psi_current < (10^(-10))] <- 10^(-10)
    psi_current[which(psi_current == max(psi_current))] <- psi_current[which(psi_current == max(psi_current))] - sum_of_psi_zero * (10^(-10))
  }

  # sample alphadp
  alphadp_proposed <- rnorm(1, alphadp_current, st_dev_for_alphadp_prop)
  while ((alphadp_proposed < 0.0001) | (alphadp_proposed > right_limit_alphadp)) {
    alphadp_proposed <- rnorm(1, alphadp_current, st_dev_for_alphadp_prop)
  }

  MH_ratio_numerator <- sum(log(psi_current^(alphadp_proposed - 1)))
  MH_ratio_denominator <- sum(log(psi_current^(alphadp_current - 1)))

  for_truncation_1 <- pnorm(right_limit_alphadp, alphadp_proposed, st_dev_for_alphadp_prop) - pnorm(0.0001, alphadp_proposed, st_dev_for_alphadp_prop)
  for_truncation_2 <- pnorm(right_limit_alphadp, alphadp_current, st_dev_for_alphadp_prop) - pnorm(0.0001, alphadp_current, st_dev_for_alphadp_prop)

  ratio_of_qproposals_num <- dnorm(alphadp_current, alphadp_proposed, st_dev_for_alphadp_prop) / for_truncation_1
  ratio_of_qproposals_den <- dnorm(alphadp_proposed, alphadp_current, st_dev_for_alphadp_prop) / for_truncation_2
  ratio_of_qproposals <- ratio_of_qproposals_num / ratio_of_qproposals_den

  MH_ratio <- exp(MH_ratio_numerator - MH_ratio_denominator) * ratio_of_qproposals
  accept_prob <- min(c(1, MH_ratio))
  if (runif(1, 0, 1) < accept_prob) {
    alphadp_current <- alphadp_proposed
  }

  # sample z
  for (i_sub in 1:N) {
    prob_for_z <- rep(0, C)
    for (i_cluster in 1:C) {
      prob_for_z[i_cluster] <- psi_current[i_cluster]
      for (i_for_prob_for_z in 1:P) {
        prob_for_z[i_cluster] <- prob_for_z[i_cluster] * phi_current[i_cluster, i_for_prob_for_z, x[i_sub, i_for_prob_for_z] + 1]
      }
    }
    prob_for_z <- prob_for_z / sum(prob_for_z)
    z_new[i_sub] <- which(as.vector(rmultinom(1, 1, prob_for_z) == 1))
  }

  z_current <- z_new
  # save samples
  if (i_iter > burnin) {
    z[i_iter - burnin, ] <- z_current
    phi[i_iter - burnin, , , ] <- phi_current
    psi[i_iter - burnin, ] <- psi_current
    alphadp[i_iter - burnin] <- alphadp_current
  }
}
