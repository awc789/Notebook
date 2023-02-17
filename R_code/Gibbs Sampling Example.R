ddirichlet = function(a,b,c) {
    x1 <- rgamma(1, a, 1)
    x2 <- rgamma(1, b, 1)
    x3 <- rgamma(1, c, 1)
    s <- x1 + x2 + x3
    d <- c(x1/s, x2/s ,x3/s)
    dim(d) = c(1,3)
    return(d)
}

sample_size <- 5000
z1 <- rep(0, sample_size)
z2 <- rep(0, sample_size)
ptheta <- matrix(0, nrow = sample_size, ncol = 3)

# z1 = rep(0, 2250)
# z2 = rep(0, 2250)
# ptheta = matrix(0, nrow=2250,ncol=3)



y1 <- 89
y2 <- 642
y3 <- 195
y4 <- 657

z1[1] <- 500
z2[1] <- 500
ptheta[1,] <- c(1/3, 1/3, 1/3)
# ptheta[1,]=c(.33,.33,.33)

a1 <- 1
a2 <- 1
a3 <- 1

for ( i in 2:sample_size) {
  # z1[i] <- rbinom(1, y2, (ptheta[i-1, 1]^2) / ((ptheta[i-1, 1]^2) + 2*ptheta[i-1, 1]*ptheta[i-1, 3]) )
  # z2[i] <- rbinom(1, y3, (ptheta[i-1, 1]^2) / ((ptheta[i-1, 1]^2) + 2*ptheta[i-1, 1]*ptheta[i-1, 3]) )
  
  z1[i] = rbinom(1, y2, (ptheta[i-1,1]^2)/  ((ptheta[i-1,1]^2)+ 2*ptheta[i-1,1]*ptheta[i-1,3])  )
  z2[i] = rbinom(1, y3, (ptheta[i-1,2]^2)/  ((ptheta[i-1,2]^2)+ 2*ptheta[i-1,2]*ptheta[i-1,3])  )
  
  m1 <- y1 + y2 + z1[i]
  m2 <- y1 + y3 + z2[i]
  m3 <- y2 + y3 - z1[i] - z2[i] + 2*y4

  ptheta[i,] <- ddirichlet(m1+a1, m2+a2, m3+a3) 
}

mean(ptheta[(sample_size/2) : sample_size, 1])