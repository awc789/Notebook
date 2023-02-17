library(MASS) # Simulate from a Multivariate Normal Distribution
library(ggplot2)
library(dirichletprocess) # Dirichlet Process
library(dbscan) # Density-based Clustering -- DBSCAN
library(dplyr)
library(class)

setwd('/Users/awc789/Library/CloudStorage/OneDrive-UniversityofStAndrews/R语言练习代码/DIB-C_data/')


name_x <- paste("x_", 1:8, sep = "")
name_workers <- paste("worker_", 1:20, sep = "")


name_workers <- paste("worker_", 1:1, sep = "")
num_dp <- 2400
MCMC_mum <- 1500

for(current_worker in name_workers){
  data <- data.frame()
  for (current_x in name_x) {
    f <- paste(current_worker, '/', current_x, '_', current_worker, '.txt', sep = "")
    data_temp <- read.table(f, sep = ',', header = TRUE)
    data <- rbind(data, data_temp)
  }
  
  # --------------------------------------------
  # DP for the clustering
  data <- data[sample(1:nrow(data)),]
  data_dp <- data[1:num_dp,]
  data_dp <- matrix(c(data_dp$X1, data_dp$X2), ncol = 2, byrow = FALSE)
  
  dpCluster <-  DirichletProcessMvnormal(data_dp)
  dpCluster <- Fit(dpCluster, MCMC_mum, progressBar = FALSE)
  plot(dpCluster)
  
  data_dp_cluster <- data[1:num_dp,]
  data_dp_cluster$clusterLabels <- dpCluster$clusterLabels
  
  # kNNdistplot(x_mix_dp, k=4)
  # abline(h=0.8, col="red")
  # db = dbscan(x_mix_dp, 0.8, 6)
  # hullplot(x_mix_dp, db$cluster)
  
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
  
  
  folder <- paste('./Center_point/', sep = "")
  # check the folder
  if(!file.exists(folder)){
    dir.create(folder)
  }
  
  f <- paste('Center_point/', 'Center_point_', current_worker, '.txt', sep = "")
  write.table(center_point, file = f, sep = ",")
  
}


