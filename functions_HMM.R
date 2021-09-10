#######LIBRARIES
library(ggpubr)
library(ggplot2)
library(Rlab)
library(HMMpa)
library(matrixcalc)

###### FUNCTION 1
###function that simulates a discrete two-state Markov Chain

#takes as input the probability that the initial chain starts at state 0,
# the probability P_00, P_11 and the number of steps

#outputs the hidden states over the time

simulate_MC <- function(s0, P_00, P_11, t){
  #load library to use rbern function
  library(Rlab)
  # set the states' names
  states = c("0", "1")
  # initial state probabilities summed to 1
  s1 = 1-s0
  
  # set the transition matrix
  # row probabilities of the matrix summed to 1
  P_01 = 1-P_00
  P_10 = 1-P_11
  P = matrix(c(P_00, P_01, P_10, P_11),
             nrow=2,byrow=TRUE,dimname=list(states,states))
  
  # create the vector to store the states at each step
  X <- vector(length = t+1)
  
  # initial state at time t=0
  # using probability of success s1
  X[1] = rbern(1, s1)
  
  #based on the current state choose the probability of success so as to find
  # the next state using Bernoulli distribution 
  for (i in 1:t){
    if (X[i] == 0){
      p = P_01
      X[i+1] = rbern(1, p)
    }
    else{
      p = P_11
      X[i+1] = rbern(1, p)
    }
  }
  return(X)
}

####### Function 2 and Function 3 produce the results without any noise added

###### FUNCTION 2
###simulates the realizations - the emitted photons in each step
#according to the state in that time-step

#takes as input the X from the simulate_MC function and the signal lambda for
# the Poisson distribution

simulate_time_series <- function(X, lambda){
  
  total = length(X)
  Y = vector(length = total)
  for (i in 1:total){
    Y[i] = rpois(1, lambda*X[i])
  }
  return(Y)
}

###### FUNCTION 3
###function that estimates the parameters P_00, P_01, P_10, P_11, s_0, s_1,
# given the vector of realizations Y
# it uses the Baum-Welch algorithm to estimate those parameters

# interested about gamma(transition matrix) and 
#delta(initial values for the marginal probability distribution of the 2 states)

#install.packages("HMMpa")
library(HMMpa)

estimate_parameters <- function(Y, s0, P_00, P_11, l){
  
  library(HMMpa)
  
  P_01 = 1 - P_00
  P_10 = 1 - P_11
  s1 = 1-s0
  
  B_W = Baum_Welch_algorithm(x=Y, m=2, delta = c(s0,s1), 
                             gamma = matrix(c(P_00,P_01,P_10,P_11), nrow=2, ncol=2),
                             distribution_class = "pois",
                             distribution_theta = list(lambda=c(0,l)))
  
  my_list = list("P_00_hat"=round(B_W$gamma[1,1],3), 
                 "P_01_hat"=round(B_W$gamma[1,2],3),
                 "P_10_hat"=round(B_W$gamma[2,1],3),
                 "P_11_hat"=round(B_W$gamma[2,2],3),
                 "s0_hat"=round(B_W$delta[1],3),
                 "s1_hat"=round(B_W$delta[2],3))
  
  return(my_list)
  
}

######Function 4
## Run the markov model many times to get different sets of simulations using
#simulate_MC function, realizations using simulate_time_series function and
#thus different estimates for the parameters using the estimate_parameters 
#function

simulate_and_estimate <- function(s0, P_00, P_11, t, lambda, iter){
  
  steps = t+1
  
  X <- matrix(NA, nrow=steps, ncol=iter)
  Y <- matrix(NA, nrow=steps, ncol=iter)
  P_00_hat = seq()
  P_01_hat = seq()
  P_10_hat = seq()
  P_11_hat = seq()
  s0_hat = seq()
  s1_hat = seq()
  
  for (i in 1:iter){
    
    X[,i] = simulate_MC(s0, P_00, P_11, t) 
    
    Y[,i] = simulate_time_series(X[,i], lambda)
    
    estimations = estimate_parameters(Y[,i], s0, P_00, P_11, lambda)
    P_00_hat[i] = estimations$P_00_hat
    P_01_hat[i] = estimations$P_01_hat
    P_10_hat[i] = estimations$P_10_hat
    P_11_hat[i] = estimations$P_11_hat
    s0_hat[i] = estimations$s0_hat
    s1_hat[i] = estimations$s1_hat
    
  }
  my_list = list("P_00_hat"=P_00_hat, "P_01_hat"=P_01_hat, "P_10_hat"=P_10_hat,
                 "P_11_hat"=P_11_hat, "s0_hat"=s0_hat, "s1_hat"=s1_hat)
  return(my_list)
}


######Function 5
# Mean Squared Error
mse <- function(actual, predicted){
  total = length(predicted)
  sum = 0
  for (i in 1:total){
    sum = sum + (actual - predicted[i])^2
  }
  return(sum/total)
}

####### Function 6, Function 7 and Function 8 produce the results with added 
#noise

###### FUNCTION 6
###function that calculates the realizations - the emitted photons in each step
#according to the state in that time-step

#takes as input the X from the simulate_MC function and the signal lambda for
# the poisson distribution
#calculates the realizations using the Poisson distribution with signal lambda 
#and noise lambda_n

simulate_time_series_noise <- function(X, lambda, lambda_n){
  
  total = length(X)
  Y = vector(length = total)
  for (i in 1:total){
    Y[i] = rpois(1, lambda*X[i]) + rpois(1, lambda_n)
  }
  return(Y)
}

###### FUNCTION 7
###Estimates the parameters P_00, P_01, P_10, P_11, s_0, s_1,
# given the vector of realizations Y with noise
# it uses the Baum-Welch algorithm to estimate those parameters

# interested about gamma(transition matrix) and 
#delta(initial values for the marginal probability distribution of the 2 states)

#install.packages("HMMpa")
library(HMMpa)

estimate_parameters_noise <- function(Y, s0, P_00, P_11, l, l_n){
  
  library(HMMpa)
  
  means = l + l_n
  P_01 = 1 - P_00
  P_10 = 1 - P_11
  s1 = 1-s0
  
  B_W = Baum_Welch_algorithm(x=Y, m=2, delta = c(s0,s1), 
                             gamma = matrix(c(P_00,P_01,P_10,P_11), nrow=2, ncol=2),
                             distribution_class = "pois",
                             distribution_theta = list(lambda=c(l_n, means)))
  
  my_list = list("P_00_hat"=round(B_W$gamma[1,1],3), 
                 "P_01_hat"=round(B_W$gamma[1,2],3),
                 "P_10_hat"=round(B_W$gamma[2,1],3),
                 "P_11_hat"=round(B_W$gamma[2,2],3),
                 "s0_hat"=round(B_W$delta[1],3),
                 "s1_hat"=round(B_W$delta[2],3),
                 "lambda_hat" = round(B_W$distribution_theta$lambda))
  
  return(my_list)
}

######Function 8
## Run the markov model many times to get different sets of simulations using
#simulate_MC function, realizations using simulate_time_series_noise function 
#and thus different estimates for the parameters using the 
#estimate_parameters_noise function

simulate_and_estimate_noise <- function(s0, P_00, P_11, t, lambda, lambda_n,
                                        iter){
  
  steps = t+1
  
  X <- matrix(NA, nrow=steps, ncol=iter)
  Y <- matrix(NA, nrow=steps, ncol=iter)
  P_00_hat = seq()
  P_01_hat = seq()
  P_10_hat = seq()
  P_11_hat = seq()
  s0_hat = seq()
  s1_hat = seq()
  
  for (i in 1:iter){
    
    X[,i] = simulate_MC(s0, P_00, P_11, t) 
    
    Y[,i] = simulate_time_series_noise(X[,i], lambda, lambda_n)
    
    estimations = estimate_parameters_noise(Y[,i], s0, P_00, P_11, lambda, lambda_n)
    P_00_hat[i] = estimations$P_00_hat
    P_01_hat[i] = estimations$P_01_hat
    P_10_hat[i] = estimations$P_10_hat
    P_11_hat[i] = estimations$P_11_hat
    s0_hat[i] = estimations$s0_hat
    s1_hat[i] = estimations$s1_hat
    
  }
  my_list = list("P_00_hat_n"=P_00_hat, "P_01_hat_n"=P_01_hat, 
                 "P_10_hat_n"=P_10_hat, "P_11_hat_n"=P_11_hat, 
                 "s0_hat_n"=s0_hat, "s1_hat_n"=s1_hat)
  return(my_list)
}

######FUNCTION 9
#takes the probabilities P_00 and P_11 as arguments and calculates
# the Frobenius norm of the transition matrix
frob_norm <- function(P_00, P_11){
  library(matrixcalc)
  P_01 = 1-P_00
  P_10 = 1-P_11
  M = matrix(c(P_00, P_01, P_10, P_11), nrow=2, ncol=2, byrow=TRUE)
  return(frobenius.norm(M))
}

######FUNCTION 10
#Function that simulates the markov chain states and estimates the realizations
# and the transition probabilities for different noise levels.
# Takes as arguments fixed transition probabilities, signal level, the number
# of steps, the number of iterations and an interval for the noise level
# For each noise level, it creates a markov chain model, which it simulates 
# for iter times.
# Outputs: 
# -mse of the estimated transition probabilities for each noise level
# compared to the real transition probabilities
# -average estimated transition probabilities for each noise level
# -average Frobenius norm for each noise level, using the matrix from the
# difference of the estimated transition matrix and the real transition matrix

investigate_noise_level <- function(s0, P_00, P_11, lambda, t, iter, 
                                    noise_min, noise_max){
  
  P_01 = 1 - P_00
  P_10 = 1 - P_11
  
  #set null vectors to store the mse between the estimated transition
  # probabilities and the true transition probabilities for each noise level
  mse_P_00_hat_level = seq()
  mse_P_01_hat_level = seq()
  mse_P_10_hat_level = seq()
  mse_P_11_hat_level = seq()
  
  #set null vectors to store the average value of each estimated transition 
  # probability for each noise level
  aver_P_00_hat_level = seq()
  aver_P_01_hat_level = seq()
  aver_P_10_hat_level = seq()
  aver_P_11_hat_level = seq()
  
  frb = seq()
  aver_frb = seq()
  
  #count the different noise levels 
  count = 1
  
  for (n in seq(noise_min, noise_max, 50)){
    
    lambda_n = n
    
    simulation = simulate_and_estimate_noise(s0, P_00, P_11, t, lambda, 
                                             lambda_n, iter) 
    
    #find the mse of the estimated transition probabilities for each noise level
    # compared to the real ones 
    mse_P_00_hat_level[count] = mse(P_00, simulation$P_00_hat_n)
    mse_P_01_hat_level[count] = mse(P_01, simulation$P_01_hat_n)
    mse_P_10_hat_level[count] = mse(P_10, simulation$P_10_hat_n)
    mse_P_11_hat_level[count] = mse(P_11, simulation$P_11_hat_n)
    
    #find the average of each estimated transition probability for each noise
    # level
    aver_P_00_hat_level[count] = mean(simulation$P_00_hat_n)
    aver_P_01_hat_level[count] = mean(simulation$P_01_hat_n)
    aver_P_10_hat_level[count] = mean(simulation$P_10_hat_n)
    aver_P_11_hat_level[count] = mean(simulation$P_11_hat_n)
    
    # get the Frobenius norm for each estimated transition matrix through the 
    # iterations
    for (i in 1:iter){
      differ = matrix(c(simulation$P_00_hat_n[i], simulation$P_01_hat_n[i], 
                           simulation$P_10_hat_n[i], simulation$P_11_hat_n[i]),
                         nrow=2, ncol=2, byrow=TRUE) - 
        matrix(c(P_00, P_01, P_10, P_11), nrow=2, ncol=2, byrow=TRUE)
      
      frb[i] = frob_norm(differ[1,1],differ[2,2])
    }
    
    # Average Frobenius norm of the error matrix for each noise level 
    aver_frb[count] = mean(frb)
    
    count = count+1
  }
  
  #Frobenius norm of the actual transition matrix 
  real_frb = frob_norm(P_00, P_11)
  
  my_list = list("mse_P_00_hat_level"=mse_P_00_hat_level,
                 "mse_P_01_hat_level"=mse_P_01_hat_level,
                 "mse_P_10_hat_level"=mse_P_10_hat_level,
                 "mse_P_11_hat_level"=mse_P_11_hat_level,
                 "aver_P_00_hat_level"=aver_P_00_hat_level,
                 "aver_P_01_hat_level"=aver_P_01_hat_level,  
                 "aver_P_10_hat_level"=aver_P_10_hat_level,  
                 "aver_P_11_hat_level"=aver_P_11_hat_level, 
                 "real_frb"=real_frb, "aver_frb"=aver_frb)
  return(my_list)
}

######FUNCTION 11
#####Estimating the numbers of frames in On and Off states, and thus the 
# transition probabilities that cause those t_on and t_off times
estimate_frames_probs <- function(t_on, t_off){
  
  #each frame has a duration of 0.05 sec.
  frames_on = t_on / 0.05
  frames_off = t_off / 0.05
  
  #the expectation of geometric distribution E=(1-p)/p
  p_10 = 1/(frames_on + 1) 
  p_01 = 1/(frames_off + 1)
  p_11 = 1-p_10
  p_00 = 1-p_01
  
  
  #use .format.zeros to remove trailing zeros
  my_list = list("p_00" = .format.zeros(p_00, zero.print=NULL),
              "p_01" = .format.zeros(p_01, zero.print=NULL),
              "p_10" = .format.zeros(p_10, zero.print=NULL),
              "p_11" = .format.zeros(p_11, zero.print=NULL))
  
  #use scipen to print the results without scientific printing with exponential
  options(scipen=1)
  
  return(my_list)
}

######FUNCTION 12
#Function that takes as arguments the t_on and t_off. It also takes as arguments
# the signal lambda, s0, t(steps) of the whole markov chain and the iterations.
# Also, the boundaries for the varying noise level.
# It computes the actual transition probabilities using the t_on and t_off 
# arguments using the "estimate_frames_probs" function. Then, using these probs
# and the rest of the arguments we will run the function investigate_noise_level.
# By doing this, we will have a fixed model with the noise level to be varying.
# For each noise level we will estimate the transition probabilities. Moreover,
# we will get mse values, average values for the estimates of the transition 
# probabilities and the Frobenius norms

# the t_on and t_off corresponds to the average time in the On and Off states
# for the slow, medium and fast photo-switching scenarios

smf_investigate_noise <- function(t_on, t_off, s0, t, iter, lambda, noise_min,
                                  noise_max){
  
  probs = estimate_frames_probs(t_on, t_off)
  P_00 = probs$p_00
  P_11 = probs$p_11
  
  est = investigate_noise_level(s0, P_00, P_11, lambda, t, iter, 
                                      noise_min, noise_max)
  
  my_list = list("P_00" = P_00, "P_11"=P_11,
                  "mse_P_00_hat_level"=est$mse_P_00_hat_level,
                   "mse_P_01_hat_level"=est$mse_P_01_hat_level,
                   "mse_P_10_hat_level"=est$mse_P_10_hat_level,
                   "mse_P_11_hat_level"=est$mse_P_11_hat_level,
                   "aver_P_00_hat_level"=est$aver_P_00_hat_level,
                   "aver_P_01_hat_level"=est$aver_P_01_hat_level,  
                   "aver_P_10_hat_level"=est$aver_P_10_hat_level,  
                   "aver_P_11_hat_level"=est$aver_P_11_hat_level, 
                   "real_frb"=est$real_frb, "aver_frb"=est$aver_frb)
  
  return(my_list)
}