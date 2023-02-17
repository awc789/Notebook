library(mvtnorm)
library(MASS)

# Create a 'true' data set (an easy one)
.create.data <- function(n)
{
  l <- list()
  l[[1]] <- list(component=1,
                 mixing.weight=0.5,
                 means=c(0,0),
                 cov=matrix(c(1,0,0,1), ncol=2, byrow=T))
  l[[2]] <- list(mixing.weight=0.3,
                 component=2,
                 means=c(5,5),
                 cov=matrix(c(1, 0.5, 0.5, 1), ncol=2, byrow=T))
  l[[3]] <- list(mixing.weight=0.2,
                 component=3,
                 means=c(10,10),
                 cov=matrix(c(1,0.75,0.75,1), ncol=2, byrow=T))
  
  do.call("rbind",sapply(l, function(e) {
    dat <- mvtnorm::rmvnorm(e$mixing.weight * n, e$means, e$cov)
    cbind(component=e$component,
          x1=dat[,1],
          x2=dat[,2])
  }))
}

# Function for covariance update
.cov <- function(n, r, dat, m, N.k)
{
  (t(r * (dat[,2:3] -m))  %*%  (( dat[,2:3]-m))) / N.k
}

# Generate starting values for means/covs/mixing weights
.init <- function()
{
  l <- list()
  l[[1]] <- list(mixing.weight=0.1,
                 means=c(-2, -2),
                 cov=matrix(c(1,0,0,1), ncol=2, byrow=T))
  l[[2]] <- list(mixing.weight=0.1,
                 means=c(10, 0),
                 cov=matrix(c(1,0,0,1), ncol=2, byrow=T))
  l[[3]] <- list(mixing.weight=0.8,
                 means=c(0, 10),
                 cov=matrix(c(1,0,0,1), ncol=2, byrow=T))
  
  l
}

# Plot the 2D contours of the estimated Gaussian components
.contour <- function(means, cov, l)
{
  X <- mvtnorm::rmvnorm(1000, means, cov)
  z <- MASS::kde2d(X[,1], X[,2], n=50)
  contour(z, drawlabels=FALSE, add=TRUE, lty=l, lwd=1.5)
}

# Do a scatter plot
.scatter <- function(dat, clusters)
{
  plot(dat[,2], dat[,3],
       xlab="X", ylab="Y", main="Three component Gaussian mixture model",
       col=c("blue", "red", "orange", "black")[clusters],
       pch=(1:4)[clusters])
  col    <- c("blue", "red", "orange")
  pch    <- 1:3
  legend <- paste("Cluster", 1:3)
  if (clusters == 4) {
    col <- "black"
    pch <- 4
    legend = "No clusters"
  }
  legend("topleft", col=col, pch=pch, legend=legend)
}

n   <- 10000
# create data with n samples
dat <- .create.data(n)
repeat
{
  # set initial parameters
  l    <- .init()
  # plot initial data
  .scatter(dat, 4)
  invisible(lapply(1:3, function(e) .contour(l[[e]]$means, l[[e]]$cov, 1)))
  # Usually we would do a convergence criterion, e.g. compare difference of likelihoods
  # but this will suffice for the hands on
  for (i in seq(50))
  {
    
    ### E step
    
    # Compute the sum of all responsibilities (for normalization)
    r <- sapply(l, function(r)
    {
      r$mixing.weight *  mvtnorm::dmvnorm(dat[,2:3], r$means, r$cov)
    })
    r <- apply(r, 1, sum)
    # Compute the responsibilities for each sample
    rs <- sapply(l, function(e)
    {
      e$mixing.weight * mvtnorm::dmvnorm(dat[,2:3], e$means, e$cov) / r
    })
    # Compute number of points per cluster
    N.k <- apply(rs, 2, sum)
    
    
    ### M step
    
    # Compute the new means
    m <- lapply(1:3, function(e)
    {
      apply(rs[,e] * dat[,2:3], 2, sum) / N.k[e]
    })
    # Compute the new covariances
    c <- lapply(1:3, function(e)
    {
      .cov(n, rs[,e], dat, m[[e]], N.k[e])
    })
    # Compute the new mixing weights
    mi <- N.k / n
    # Update the old parameters
    l <- lapply(1:3, function(e)
    {
      list(mixing.weight = mi[e], means=m[[e]], cov=c[[e]])
    })
    # Plot a 2D density (contour) to show the estimated means and covariances
    if (i %% 5 == 0)
    {
      Sys.sleep(1.5)
      .scatter(dat, apply(rs, 1, which.max))
      invisible(lapply(1:3, function(e) .contour(l[[e]]$means, l[[e]]$cov, e + 1)))
    }
  }
}