library(parallel)
library(MASS)
library(dirichletprocess) # Dirichlet Process
library(doParallel)
library(foreach)


# -------------------------------------------------------------------------
path <- paste('/Users/awc789/Library/CloudStorage/OneDrive-UniversityofStAndrews/GitHub/Notebook/R_code', '/DC_data/', sep = '')

if(!file.exists(path)){
  dir.create(path)
}
setwd(path)

name_x <- paste("x_", 1:8, sep = "")
name_GMM <- paste("GMM_", 1:4, sep = "")
workers_num <- detectCores() * 2
name_workers <- paste("worker_", 1:workers_num, sep = "")

num_dp <- 100 # 放入DP的数据量
MCMC_mum <- 10 # 迭代次数

parallel_function <- function(index){
  # 'original' or 'GMM'
  p <- 'original'
  
  dp_num <- 3200    # The number of subsamping for DP clustering
  MCMC_num <- 1500  # The number of samples obtained during MCMC
  
  current_worker <- paste("worker_", index, sep = "")
  name_GMM <- paste("GMM_", 1:4, sep = "")
  name_x <- paste("x_", 1:8, sep = "")
  
  
  data <- data.frame()
  if (p == 'original') {
    for (current_x in name_x) {
      f <- paste(path, current_worker, '/', current_x, '_', current_worker, '.txt', sep = "")
      data_temp <- read.table(f, sep = ',', header = TRUE)
      data <- rbind(data, data_temp)
    }
  } else if (p == 'GMM') {
    for (current_GMM in name_GMM) {
      f <- paste(path, current_worker, '/', current_GMM, '_', current_worker, '.txt', sep = "")
      data_temp <- read.table(f, sep = ',', header = TRUE)
      data <- rbind(data, data_temp)
    }
  }
  
  # --------------------------------------------
  # DP for the clustering
  data <- data[sample(1:nrow(data)),] # shuffle
  data_dp <- data[1:dp_num,]
  data_dp <- matrix(c(data_dp$X1, data_dp$X2), ncol = 2, byrow = FALSE)
  
  
  dpCluster <-  DirichletProcessMvnormal(data_dp)
  dpCluster <- Fit(dpCluster, MCMC_num, progressBar = FALSE)
  #plot(dpCluster)
  
  data_dp_cluster <- data[1:dp_num,]
  data_dp_cluster$clusterLabels <- dpCluster$clusterLabels
  
  label_result <- t(data.frame(table(dpCluster$clusterLabels)))
  drop_clusterlabels <- c()
  for (i in seq(ncol(label_result))) {
    if(label_result[2,i] <= (nrow(data_dp_cluster)/ncol(label_result))){
      drop_clusterlabels <- append(drop_clusterlabels, i)
      data_dp_cluster <- data_dp_cluster[data_dp_cluster$clusterLabels != i,]
    }
  }
  
  remain_clusterlabels <- data.frame(table(data_dp_cluster$clusterLabels))
  center_point <- data.frame(X1 = c(1:nrow(remain_clusterlabels)), 
                             X2 = c(1:nrow(remain_clusterlabels)), 
                             clusterlabels = c(1:nrow(remain_clusterlabels)))
  
  remain_clusterlabels <- remain_clusterlabels$Var1
  for(i in remain_clusterlabels){
    i <- as.numeric(i)
    temp_X1 <- data_dp_cluster[data_dp_cluster$clusterLabels == i,]$X1
    temp_X2 <- data_dp_cluster[data_dp_cluster$clusterLabels == i,]$X2
    
    center_point[i,1] <- mean(temp_X1)
    center_point[i,2] <- mean(temp_X2)
    center_point[i,3] <- i
  }
  
  folder <- paste(path, '/Center_point/', sep = "")
  # check the folder
  if(!file.exists(folder)){
    dir.create(folder)
  }
  
  f <- paste(path, 'Center_point/Center_point_', current_worker, '.txt', sep = "")
  write.table(center_point, file = f, sep = ",")
  
  return('Null')
}

# -------------------------------------------------------------------------

# cl_cores <- detectCores(logical = F)
# cl <- makeCluster(cl_cores)
# registerDoParallel(cl)
# clusterExport(cl,deparse(substitute(parallel_function), DirichletProcessMvnormal))

# x <- foreach(i = 1:workers_num, .combine <- rbind, .packages = c('dirichletprocess')) %dopar% parallel_function(i)

# stopCluster(cl)

# -------------------------------------------------------------------------

time_1 <- function(){
  cl_cores <- detectCores(logical = F)
  cl <- makeCluster(cl_cores)
  registerDoParallel(cl) 
  
  func <- function(ii){
    Sys.sleep(1)
    #print(ii)
  }
  x <- foreach(ii=1:16) %dopar% func(ii)
  
  stopCluster(cl)
}

time_2 <- function(){
  func <- function(ii){
    Sys.sleep(1)
  }
  ii <- c(1:16)
  for (i in ii) {
    func(ii)
  }
}

system.time(time_1())
system.time(time_2())
