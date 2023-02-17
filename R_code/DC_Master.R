library(ggplot2)
library(dirichletprocess) # Dirichlet Process
library(dbscan) # Density-based Clustering -- DBSCAN


setwd('/Users/awc789/Library/CloudStorage/OneDrive-UniversityofStAndrews/R语言练习代码/DIB-C_data/')

name_x <- paste("x_", 1:8, sep = "")
name_workers <- paste("worker_", 1:20, sep = "")

data <- data.frame()
for(current_worker in name_workers){
  f <- paste('Center_point/', 'Center_point_', current_worker, '.txt', sep = "")
  data_temp <- read.table(f, sep = ',', header = TRUE)
  data <- rbind(data, data_temp)
}

ggplot(data, aes(x = X1, y = X2, colour = clusterlabels)) + 
  geom_point(size=3)

data_central <- matrix(c(data$X1, data$X2), ncol = 2, byrow = FALSE)
# dpCluster <-  DirichletProcessMvnormal(data_dp)
# dpCluster <- Fit(dpCluster, 2000, progressBar = FALSE)
# plot(dpCluster)


# kNNdistplot(data_dp, k=4)
# abline(h=0.3, col="red")
db = dbscan(data_central, 5, 4)
hullplot(data_central, db$cluster)

data$clusterlabels <- db$cluster
ggplot(data, aes(x = X1, y = X2, colour = clusterlabels)) + 
  geom_point(size=3)

center_point <- data.frame(X1 = c(1:nrow(data.frame(table(data$clusterlabels)))), 
                           X2 = c(1:nrow(data.frame(table(data$clusterlabels)))), 
                           clusterlabels = c(1:nrow(data.frame(table(data$clusterlabels)))))
for(i in data.frame(table(data$clusterlabels))$Var1){
  i <- as.numeric(i)
  temp_X1 <- data[data$clusterlabels == i,]$X1
  temp_X2 <- data[data$clusterlabels == i,]$X2
  
  center_point[i,1] <- mean(temp_X1)
  center_point[i,2] <- mean(temp_X2)
  center_point[i,3] <- i
}

master_clustering <- function(data, center){
  data_result <- data.frame()
  for(i in seq(nrow(data))){
    X1 <- data[i,]$X1
    X2 <- data[i,]$X2
    data_result[i,1] <- X1
    data_result[i,2] <- X2
    for(j in seq(nrow(center))){
      distance <- sqrt((X1 - center[j,]$X1)^2 + (X2 - center[j,]$X2)^2)
      data_result[i,j+2] <- distance
    }
  }
  k <- ncol(data_result)
  # distace_name <- paste("", 1:(k-2), sep = "")
  distace_name <- c(1:(k-2))
  colnames(data_result) <- c('X1', 'X2', distace_name)
  data_result$clusterlabels <- colnames(data_result[,3:k])[apply(data_result[,3:k],1,which.min)]
  return(data_result)
}

data_master <- data.frame()
for(current_x in name_x){
  f <- paste(current_x, '_total.txt', sep = "")
  data_temp <- read.table(f, sep = ',', header = TRUE)
  data_master <- rbind(data_master, data_temp)
}

sample_size <- 8000
data_master <- data_master[sample(1:nrow(data_master)), ]

data_master_sample <- data_master[1:sample_size, ]

master_result <- master_clustering(data_master_sample, center_point)

ggplot(master_result, aes(x = X1, y = X2, colour = clusterlabels)) + 
  geom_point(size=3)


