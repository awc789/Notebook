library(parallel)
library(MASS)
library(dirichletprocess) # Dirichlet Process


# -------------------------------------------------------------------------
path <- paste('/Users/awc789/Library/CloudStorage/OneDrive-UniversityofStAndrews/R语言练习代码', '/DC_data_v2/', sep = '')

if(!file.exists(path)){
  dir.create(path)
}
setwd(path)

name_workers <- paste("worker_", 1:1, sep = "")
name_GMM <- paste("GMM_", 1:4, sep = "")

num_dp <- 100 # 放入DP的数据量
MCMC_mum <- 10 # 迭代次数

current_worker <- name_workers[1]
data <- data.frame()

# 读取数据
for (current_GMM in name_GMM) {
  f <- paste(path, current_worker, '/', current_GMM, '_', current_worker, '.txt', sep = "")
  data_temp <- read.table(f, sep = ',', header = TRUE)
  data <- rbind(data, data_temp)
}

data <- data[sample(1:nrow(data)),] # 打乱顺序
data_dp <- data[1:num_dp,] # 获得DP数据
data_dp <- matrix(c(data_dp$X1, data_dp$X2), ncol = 2, byrow = FALSE)

dpCluster <-  DirichletProcessMvnormal(data_dp)
dpCluster <- Fit(dpCluster, MCMC_mum, progressBar = FALSE)


# Parallel Operation ------------------------------------------------------

r <- mclapply(1:10, function(i) {Sys.sleep(10)}, mc.cores = 5) # 20s
r <- mclapply(1:10, function(i) {Sys.sleep(10)}, mc.cores = 1) # 100s

fx <- function(dpCluster) {
  result <- Fit(dpCluster, MCMC_mum, progressBar = FALSE)
  return(result)
  }

dpCluster <-  DirichletProcessMvnormal(data_dp)
system.time({
  dpCluster <- Fit(dpCluster, MCMC_mum, progressBar = FALSE)
  })


dpCluster <-  DirichletProcessMvnormal(data_dp)
system.time({
  dpCluster <- mclapply(dpCluster, fx(dpCluster), mc.cores = 6)
  
})


plot(dpCluster)

# -------------------------------------------------------------------------


