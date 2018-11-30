rm(list = ls())

setwd('E:\\Fall Term\\Data Mining - MGMT 57100\\Project\\A Balanced kMeans Algorithm\\D3 Wrestling Deliverables')
getwd()

library(lpSolve)
library(dummies)


school_coordinates <- read.csv("School_Coordinates.csv")
Power_Ratings <- read.csv("Power_Ratings.csv")
school_coordinates <- merge(school_coordinates, Power_Ratings, by.x = 'ï..School', by.y = 'School')

data <- school_coordinates[, c(2:4)]
data$NTTP_Score[is.na(data$NTTP_Score)] <- mean(data$NTTP_Score, na.rm = T)
k <- 3

KMeans_default_v2 <- function (data = data, k = k, standardize = FALSE, max_ele = NA, min_ele = NA, max_iter = 20, seed = 42, target_var = NA, delta = 1) {
  
  data <- data.frame(data)
  
  # Standardize (if neccessary)
  if (standardize) {
    for (i in 1:ncol(data)) {
      #i<-1
      sigma_val <- sd(as.numeric(as.vector(data[,i])))
      mean_val <- mean(as.numeric(as.vector(data[,i])))
      if (sigma_val != 0) {
        data[,colnames(data)[i]] <- (as.numeric(as.vector(data[,i])) - mean_val)/sigma_val
      } else {
        warning('The variable does not posses any variation to standardize')
      }
    }
  }
  
  if (!is.na(target_var)) {
    target_data <- data[, target_var]
    data[, target_var] <- NULL
  }
  
  if(ncol(data) < 2) {
    stop('Atleast provide two or more columns for clustering. For 1-D use cut funtion of base R')
  }
  
  center_index <- sample(1:nrow(data), k)
  centers <- data[center_index,!(colnames(data) %in% c('cluster'))]
  
  data$cluster <- 0
  iter <- 1
  balanced_flag <- FALSE
  cluster_list <- data.frame(c(1:nrow(data)))
  
  while (iter < max_iter) {
    print(paste('iteration:', iter))
    cluster_dist <- vector(mode = 'double', length = k)
    cluster_dist_matrix <- NULL
    for (i in 1:nrow(data)) {
      for (j in 1:k) {
        cluster_dist[j] <- dist(rbind(data[i,!(colnames(data) %in% c('cluster'))], centers[j,!(colnames(data) %in% c('cluster'))]))
      }
      cluster_dist_matrix <- cbind(cluster_dist_matrix, cluster_dist)
      data[i, 'cluster'] <- which.min(cluster_dist)
    }
    
    
    if (balanced_flag) {
      #balanced_flag <- FALSE
      print('balancing clusters')
      
      na <- data.frame(table(data$cluster))
      na <- na$Freq
      
      if (!min_ele | !is.numeric(min_ele)) {
        min_ele <- floor(nrow(data)/k)
      }
      
      
      if (!max_ele | !is.numeric(max_ele)) {
        max_ele <- ceiling(nrow(data)/(k-1))
      }
      
      
      if (min(na) >= min_ele & max(na) <= max_ele) {
        print('clusters balanced')
        iter <- max_iter
      } else {
        #print(na)
        #print(c(min(na), max(na), min_ele, max_ele))
        b <- matrix(0, nrow = nrow(data), ncol = k)
        
        cluster_dist_matrix <- cluster_dist_matrix/max(cluster_dist_matrix)
        
        objective.in  <- matrix((-1*(t(cluster_dist_matrix))), ncol = 1)
        
        col_cond <- matrix(diag(k), k, k*nrow(data))
        col_choice <- NULL
        for (i in 1:k) {
          col_choice <- c(col_choice, rep(i,nrow(data)))
        }
        col_cond <- diag(k)
        col_cond <- col_cond[,col_choice]
        
        const.mat <- rbind(diag(nrow = k*nrow(data)), diag(nrow = k*nrow(data)), matrix(diag(nrow(data)), nrow(data), k*nrow(data)), col_cond, col_cond)
        #nrow(const.mat)
        #ncol(const.mat)
        
        const.dir  <- c(rep('<=', k*nrow(data)), rep('>=', k*nrow(data)), rep('=', nrow(data)), rep('>=', k), rep('<=', k))
        
        const.rhs <- c(rep(1, k*nrow(data)), rep(0, k*nrow(data)), rep(1, nrow(data)), rep(min_ele, k), rep(max_ele, k))
        
        if(!is.na(target_var)) {
          
          avg_ele <- mean(min_ele, max_ele)
          min_const <- avg_ele*mean(target_data) - delta*sd(target_data)
          max_const <- avg_ele*mean(target_data) + delta*sd(target_data)
          
          s_col_cond <- sweep(col_cond, MARGIN = 2, target_data, '*')
          
          const.mat <- rbind(const.mat, s_col_cond, s_col_cond)
          const.dir <- c(const.dir, rep('>=', k), rep('<=', k))
          const.rhs <- c(const.rhs, rep(min_const, k), rep(max_const, k))
          
        } else {
          min_const <- NA
          max_const <- NA
        }
        
        optimum <-  lp(direction="max",  objective.in, const.mat, const.dir,  const.rhs)
        
        if(optimum$status) {
          warnings("Couldn't find a optimized balaced class")
        }
        
        best_sol <- optimum$solution
        best_sol <- matrix(best_sol, ncol = k)
        data$cluster <- max.col(best_sol)
        print(c(min_ele, max_ele, min_const, max_const))
        
      }
    }
    
    center_past <- centers
    for (j in 1:k) {
      centers[k,] <- colMeans(data[data$cluster == k,!(colnames(data) %in% c('cluster'))])
    }
    
    if (iter == max_iter) {
      warning('Reached last iteration')
    }
    if (identical(centers, center_past)) {
      print('clusters generated')
      if (balanced_flag) {
        print('clusters generated and balanced')
        iter <- max_iter
      } else {
        print('now balancing the clusters')
        balanced_flag = TRUE
        iter <- iter + 1
      }
    } else {
      iter <- iter + 1
      print('still looping')
      plot(data$Longitude, data$Latitude,col = rainbow(k)[data$cluster])
      
    }
  }
  cluster_list <- cbind(cluster_list, data$cluster)
  #write.csv(cluster_list, 'cluster_list.csv')
  return(list(cluster = data$cluster, centers = centers))
}
k <- 6
#clustered_data <- KMeans_default_v2(data, k, seed = 53)
clustered_data <- KMeans_default_v2(data, k, seed = 58, target_var = 'NTTP_Score', delta = 1.5)
plot(data$Longitude, data$Latitude, col = rainbow(k)[clustered_data$cluster], xlab = 'Latitude', ylab = 'Longitude', main = paste('Clustering result for k =', k))
library(dplyr)

data$cluster <- clustered_data$cluster
write.csv(data, 'univ_clusters.csv')

trad_kmeans <- kmeans(x=data[, c(1:2)], centers=k, iter.max = 20, nstart = 1) 
data$trad_cluster <- trad_kmeans$cluster
write.csv(data, 'univ_clusters_with_trad.csv')

data %>% mutate(cluster = clustered_data$cluster) %>% group_by(cluster) %>% summarise(mean_power = mean(NTTP_Score))
data %>% mutate(trad_cluster = data$trad_cluster) %>% group_by(trad_cluster) %>% summarise(mean_power = mean(NTTP_Score))
