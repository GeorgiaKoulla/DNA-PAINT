library(ggpubr)
library(ggplot2)
library(Rlab)
library(HMMpa)
library(matrixcalc)

###### FUNCTION 1
###function that simulates the states of a fHMM

simulate_fHMM <- function(s0, P_00, P_11, t){
  library(Rlab)

  s1 = 1-s0
  P_01 = 1-P_00
  P_10 = 1-P_11
  
  K = length(s0)
  X = matrix(NA, nrow=K, ncol = t+1)
  
  for (i in 1:K){
    X[i,1] <- rbern(1, s1[i])
    
    for (j in 2:(t+1)){
      if (X[i,1] == 0){
        p = P_01[i]
        X[i, j] = rbern(1, p)
      }
      else{
        p = P_11[i]
        X[i, j] = rbern(1, p)
      }
    }
  }
  return(X)
}

###### FUNCTION 2
###function that calculates the realizations(without noise)
# we can add the states and then use the Poisson distribution, because the states
# are independent.

simulate_time_series_fHMM <- function(X, lambda){ 
  
  K = length(lambda)
  total = ncol(X)
  Y = vector(length = total)
  for (j in 1:total){
    sum = 0
    for (i in 1:K){
     sum = sum + rpois(1, lambda[i]*X[i,j]) 
    }
    Y[j] = sum
  }
  return(Y)
}

###### FUNCTION 3
###function that calculates the realizations with noise
simulate_time_series_noise_fHMM <- function(X, lambda, lambda_n){ #give C (covariance) and means
  
  K = length(lambda)
  total = ncol(X)
  Y = vector(length = total)
  for (j in 1:total){
    sum = 0
    for (i in 1:K){
      sum = sum + rpois(1, lambda[i]*X[i,j]) + rpois(1, lambda_n)
    }
    Y[j] = sum
  }
  return(Y)
}


