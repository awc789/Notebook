library(MASS) # Simulate from a Multivariate Normal Distribution
library(ggplot2)
library(dirichletprocess) # Dirichlet Process
library(dbscan) # Density-based Clustering -- DBSCAN

## The Synthetic Data Setup
workers_num <- 20
partition <- 50000
num <- partition * workers_num

miu_matirx <- matrix(c(6, 4, 8, 22.5, 20, 22, 22, 6.5, 1.5, 6, 6, 1.5, 8, 31, 21, 29),
                     nrow = 2, ncol = 8, byrow=FALSE)

covar_matrix <-array(c(4.84, 0, 0, 2.89, 3.61, 5.05, 5.05, 14.44, 3.61, -5.05,
                       -5.05, 14.44, 12.25, 0, 0, 3.24, 3.24, 0, 0, 12.25, 
                       14.44, 0, 0, 2.25, 2.25, 0, 0, 17.64, 2.25, 4.2, 4.2, 16),
                     dim = c(2, 2, 8))

weight_GMM <- c(1/4, 1/4, 1/4, 1/4)
w1 <- c(1/3, 1/3 ,1/3)
w2 <- c(1/2, 1/2)
w3 <- c(1/2, 1/2)
w4 <- c(1)

## Generate the Initial Data
x_1 <- data.frame(mvrnorm(n = num, miu_matirx[,1], covar_matrix[, , 1]))
x_6 <- data.frame(mvrnorm(n = (num*1.5), miu_matirx[,2], covar_matrix[, , 2]))
x_4 <- data.frame(mvrnorm(n = (num*1.5), miu_matirx[,3], covar_matrix[, , 3]))
x_8 <- data.frame(mvrnorm(n = (num*3), miu_matirx[,4], covar_matrix[, , 4]))
x_2 <- data.frame(mvrnorm(n = num, miu_matirx[,5], covar_matrix[, , 5]))
x_3 <- data.frame(mvrnorm(n = num, miu_matirx[,6], covar_matrix[, , 6]))
x_7 <- data.frame(mvrnorm(n = (num*1.5), miu_matirx[,7], covar_matrix[, , 7]))
x_5 <- data.frame(mvrnorm(n = (num*1.5), miu_matirx[,8], covar_matrix[, , 8]))

## Give Data lables
x_1$component = rep(1,num)
x_2$component = rep(1,num)
x_3$component = rep(1,num)
x_4$component = rep(4,(num*1.5))
x_5$component = rep(4,(num*1.5))
x_6$component = rep(8,(num*1.5))
x_7$component = rep(8,(num*1.5))
x_8$component = rep(12,(num*3))

## Give Data sub_lables
x_1$subcomponent = rep(1,num)
x_2$subcomponent = rep(2,num)
x_3$subcomponent = rep(3,num)
x_4$subcomponent = rep(4,(num*1.5))
x_5$subcomponent = rep(5,(num*1.5))
x_6$subcomponent = rep(6,(num*1.5))
x_7$subcomponent = rep(7,(num*1.5))
x_8$subcomponent = rep(8,(num*3))

## Save the data
write.table(x_1, file = "x_1_total.txt", sep = ",")
write.table(x_2, file = "x_2_total.txt", sep = ",")
write.table(x_3, file = "x_3_total.txt", sep = ",")
write.table(x_4, file = "x_4_total.txt", sep = ",")
write.table(x_5, file = "x_5_total.txt", sep = ",")
write.table(x_6, file = "x_6_total.txt", sep = ",")
write.table(x_7, file = "x_7_total.txt", sep = ",")
write.table(x_8, file = "x_8_total.txt", sep = ",")




# ------------------------------------------------------------------------------
# Functions  v1 ----------------------------------------------------------------
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Functions  v2 ----------------------------------------------------------------
# ------------------------------------------------------------------------------

central_distance <- function(center){
  data_result <- data.frame()
  count <- 0
  for(i in seq(nrow(center))){
    X1 <- center[i,]$X1
    X2 <- center[i,]$X2
    for(j in c(i:nrow(center))){
      if(i != j){
        count <- count + 1
        distance <- sqrt((X1 - center[j,]$X1)^2 + (X2 - center[j,]$X2)^2)
        data_result[count,1] <- distance
      }
    }
  }
  colnames(data_result) <- c('distacen')
  return(min(data_result$distacen))
}

# ------------------------------------------------------------------------------

master_clustering_v2 <- function(data, center){
  min_central_distance <- central_distance(center_point) / nrow(center_point)
  
  data_result <- data.frame(t(data.frame(c(seq(nrow(center)+8)))), row.names = NULL)
  distace_name <- c(1:nrow(center))
  colnames(data_result) <- c('X1', 'X2', distace_name, 'clusterlabels', 
                             'close_to_central', 'large_density', 'Edge_point', 
                             'Edge_diff', 'final_lable')
  
  for(i in seq(nrow(data))){
    X1 <- data[i,]$X1
    X2 <- data[i,]$X2
    data_result[i,1] <- X1
    data_result[i,2] <- X2
    for(j in seq(nrow(center))){
      distance <- sqrt((X1 - center[j,]$X1)^2 + (X2 - center[j,]$X2)^2)
      data_result[i,j+2] <- distance
    }
  
    min_distance <- min(data_result[i,3:(nrow(center)+2)])
    temp <- data_result[i,3:(nrow(center)+2)] == min_distance
    label <- col(data_result[i,3:(nrow(center)+2)])[temp]
    data_result[i,(nrow(center)+3)] <- label
    
    if(min_distance <= min_central_distance){
      data_result[i,(nrow(center)+4)] <- label
    }else{
      data_result[i,(nrow(center)+4)] <- 0
    }
    data_result[i,(nrow(center)+5)] <- 0
    data_result[i,(nrow(center)+6)] <- 0
    data_result[i,(nrow(center)+7)] <- 0
    data_result[i,(nrow(center)+8)] <- 0
  }
  return(data_result)
}


# ------------------------------------------------------------------------------
# Process ----------------------------------------------------------------------
# ------------------------------------------------------------------------------
path <- paste(getwd(), '/DC_data/', sep = '')
setwd(path)

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


data_master <- data.frame()
for(current_x in name_x){
  f <- paste(current_x, '_total.txt', sep = "")
  data_temp <- read.table(f, sep = ',', header = TRUE)
  data_master <- rbind(data_master, data_temp)
}

sample_size <- 8000
data_master <- data_master[sample(1:nrow(data_master)), ]
data_master_sample <- data_master[1:sample_size, ]
# master_result <- master_clustering(data_master_sample, center_point)
master_result <- master_clustering_v2(data_master_sample, center_point)


# -------------------------------------------------------------------------
master_result$large_density <- NA
# master_result$large_density <- 0

range_value <- 0.5
range_num <- 10
for(i in seq(nrow(master_result))){
  temp <- master_result[i,]
  within_range <- ((temp$X1-range_value) < master_result$X1 & master_result$X1 < (temp$X1+range_value)) & ((temp$X2-range_value) < master_result$X2 & master_result$X2 < (temp$X2+range_value))
  if(sum(within_range)>=range_num){
    range_data <- master_result[within_range,]
    range_clusterlabels <- table(range_data$clusterlabels)
    range_close_to_central <- table(range_data$close_to_central)
    if(length(range_clusterlabels)==1 & length(range_close_to_central)==1){
      master_result[i,(nrow(center_point)+5)] <- temp$clusterlabels
    }
  }
}

sub_master_result <- master_result[!is.na(master_result$large_density),]
sub_master_result <- sub_master_result[sub_master_result$close_to_central==0,]

ggplot(sub_master_result, aes(x = X1, y = X2, colour = large_density)) + 
  geom_point(size=3)


# -------------------------------------------------------------------------
# master_result$Edge_point <- NA
master_result$Edge_point <- 0
master_result$Edge_diff <- NA
# master_result$Edge_diff <- 0

min_central_distance <- central_distance(center_point) / (nrow(center_point) * 0.75)
for(i in seq(nrow(master_result))){
  temp <- master_result[i, 3:(2+nrow(center_point))]
  order_min <- order(temp)[1:2]
  if(temp[order_min[2]] - temp[order_min[1]] <= min_central_distance){
    master_result[i,(nrow(center_point)+6)] <- (-1) * order_min[1] * order_min[2]
    master_result[i,(nrow(center_point)+7)] <- (temp[order_min[2]] - temp[order_min[1]])
  }
}

edge_master_result <- master_result[master_result$Edge_point!=0,]


ggplot(edge_master_result, aes(x = X1, y = X2, colour = Edge_point)) + 
  geom_point(size=3)


# -------------------------------------------------------------------------

data_result <- data.frame()
for(i in seq(nrow(edge_master_result))){
  temp <- edge_master_result[i,]
  X1 <- temp$X1
  X2 <- temp$X2
  for(j in seq(nrow(sub_master_result))){
    new_sub <- sub_master_result[j,]
    if(row.names(temp) == row.names(new_sub)){
      data_result[i,j] <- Inf
    }else{
      distance <- sqrt((X1 - new_sub$X1)^2 + (X2 - new_sub$X2)^2)
      data_result[i,j] <- distance
    }
  }
}

colnames(data_result) <- c(1:(nrow(sub_master_result)))
rownames(data_result) <- row.names(edge_master_result)

master_result$final_lable <- master_result$close_to_central
closest_num <- 20
percentage <- 0.75
edge_percent <- closest_num * percentage
for(i in seq(nrow(data_result))){
  row_final <- as.numeric(row.names(data_result[i,]))
  temp <- order(data_result[i,])[1:closest_num]
  label_list <- c()
  for(j in seq(closest_num)){
    label <- sub_master_result[temp[j],]$large_density
    label_list <- append(label_list, label)
  }
  if(table(label_list)[order(table(label_list))[length(table(label_list))]] >= edge_percent){
    final_label <- as.numeric(row.names(data.frame(table(label_list)[order(table(label_list))[length(table(label_list))]])))
    master_result[row_final, (nrow(center_point)+8)] <- final_label
  }else{
    master_result[row_final, (nrow(center_point)+8)] <- Inf 
  }
}

# -------------------------------------------------------------------------

for(i in seq(nrow(master_result))){
  if(master_result[i,]$final_lable == 0){
    master_result[i,]$final_lable <- master_result[i,]$clusterlabels
  }else{
    if(master_result[i,]$final_lable == Inf){
      label_1 <- master_result[i,]$clusterlabels
      label_2 <- master_result[i,]$Edge_point / (-1 * label_1)
      p <- runif (1)
      if(p <= 0.5){
        master_result[i,]$final_lable <- label_1
      }else{
        master_result[i,]$final_lable <- label_2
      }
    }
  }
}

ggplot(master_result, aes(x = X1, y = X2, colour = final_lable)) + 
  geom_point(size=3)



# -------------------------------------------------------------------------


ggplot(master_result, aes(x = X1, y = X2, colour = close_to_central)) + 
  geom_point(size=3)

ggplot(master_result, aes(x = X1, y = X2, colour = clusterlabels)) + 
  geom_point(size=3)








