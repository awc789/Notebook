library(MASS) # Simulate from a Multivariate Normal Distribution
library(ggplot2)

setwd('/Users/awc789/Library/CloudStorage/OneDrive-UniversityofStAndrews/R语言练习代码/DIB-C_data/')

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
x_6 <- data.frame(mvrnorm(n = num, miu_matirx[,2], covar_matrix[, , 2]))
x_4 <- data.frame(mvrnorm(n = num, miu_matirx[,3], covar_matrix[, , 3]))
x_8 <- data.frame(mvrnorm(n = num, miu_matirx[,4], covar_matrix[, , 4]))
x_2 <- data.frame(mvrnorm(n = num, miu_matirx[,5], covar_matrix[, , 5]))
x_3 <- data.frame(mvrnorm(n = num, miu_matirx[,6], covar_matrix[, , 6]))
x_7 <- data.frame(mvrnorm(n = num, miu_matirx[,7], covar_matrix[, , 7]))
x_5 <- data.frame(mvrnorm(n = num, miu_matirx[,8], covar_matrix[, , 8]))

GMM_1 <- w1[1] * x_1 + w1[2] * x_2 + w1[3] * x_3
GMM_2 <- w2[1] * x_4 + w2[2] * x_5
GMM_3 <- w3[1] * x_6 + w3[2] * x_7
GMM_4 <- w4[1] * x_8

Mixture_GMM <- weight_GMM[1] * GMM_1 + weight_GMM[2] * GMM_2 + weight_GMM[3] * GMM_3 + weight_GMM[4] * GMM_4


## Give Data lables
x_1$component = rep(1.1,num)
x_2$component = rep(1.2,num)
x_3$component = rep(1.3,num)
x_4$component = rep(3.1,num)
x_5$component = rep(3.2,num)
x_6$component = rep(6.1,num)
x_7$component = rep(6.2,num)
x_8$component = rep(12.1,num)


num_sample <- 200
x_mix <- rbind(x_1[1:num_sample,], x_2[1:num_sample,])
x_mix <- rbind(x_mix, x_3[1:num_sample,]) 
x_mix <- rbind(x_mix, x_4[1:num_sample,]) 
x_mix <- rbind(x_mix, x_5[1:num_sample,]) 
x_mix <- rbind(x_mix, x_6[1:num_sample,]) 
x_mix <- rbind(x_mix, x_7[1:num_sample,]) 
x_mix <- rbind(x_mix, x_8[1:num_sample,]) 

ggplot(x_mix, aes(x = X1, y = X2, colour = component)) + 
  geom_point(size=3)


# --------------------------------------------
Mixture_GMM$component = rep(-10,num)
x_mix_with_GMM <- rbind(x_mix, Mixture_GMM[1:num_sample,]) 

ggplot(x_mix_with_GMM, aes(x = X1, y = X2, colour = component)) + 
  geom_point(size=3)


write.table(x_1, file = "x_1_total.txt", sep = ",")
write.table(x_2, file = "x_2_total.txt", sep = ",")
write.table(x_3, file = "x_3_total.txt", sep = ",")
write.table(x_4, file = "x_4_total.txt", sep = ",")
write.table(x_5, file = "x_5_total.txt", sep = ",")
write.table(x_6, file = "x_6_total.txt", sep = ",")
write.table(x_7, file = "x_7_total.txt", sep = ",")
write.table(x_8, file = "x_8_total.txt", sep = ",")

name_x <- paste("x_", 1:8, sep = "")
name_workers <- paste("worker_", 1:20, sep = "")
count <- 0
for(i in list(x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8)){
  count <- count +1
  current_x <- name_x[count]
  # For the workers
  for(k in rep(1:workers_num)){
    current_worker <- name_workers[k]
    folder <- paste("/Users/awc789/Library/CloudStorage/OneDrive-UniversityofStAndrews/R语言练习代码/DIB-C_data", '/', current_worker, sep = "")
    # check the folder
    if(!file.exists(folder)){
      dir.create(folder)
    }
    # for each workers
    start_loop <- partition * (k-1) + 1
    end_loop <- partition * k
    temp <- i[start_loop:end_loop,]
    f <- paste(current_worker, '/', current_x, '_', current_worker, '.txt', sep = "")
    write.table(temp, file = f, sep = ",")
  }
}

write.table(Mixture_GMM, file = "Mixture_GMM_total.txt", sep = ",")

for(k in rep(1:workers_num)){
  current_worker <- name_workers[k]
  folder <- paste("/Users/awc789/Library/CloudStorage/OneDrive-UniversityofStAndrews/R语言练习代码/DIB-C_data", '/', current_worker, sep = "")
  # for each workers
  start_loop <- partition * (k-1) + 1
  end_loop <- partition * k
  temp <- i[start_loop:end_loop,]
  f <- paste(current_worker, '/', 'Mixture_GMM_', current_worker, '.txt', sep = "")
  write.table(temp, file = f, sep = ",")
}
