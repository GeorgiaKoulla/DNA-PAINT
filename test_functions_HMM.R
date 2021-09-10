#### USE FUNCTIONS WITHOUT NOISE

# get the simulations of the Markov chain 
# simulation contains a list of:
# P_00, P_01, P_10, P_11, s0, s1, X, t
  X = simulate_MC(0.2, 0.3, 0.8, 50)

# plot the states at each time
  library(ggplot2)
  time = seq(0, 50)
  data_x = cbind.data.frame(time, X)
  
  ggplot(data_x, aes(time, X)) + ylim(0,3) + geom_step() +
    ggtitle("Hidden States for each time-step") + xlab("Time") + 
    ylab("State")

# get the realizations of the Markov model
# time series contains a list of:
# lambda(l) and realizations Y
  Y = simulate_time_series(X, lambda = 600)

# plot the realization of photons in each step
  data_y = cbind.data.frame(time, Y)
  
  ggplot(data_y, aes(time, Y)) + geom_step() +
    ggtitle("Photo-emissions in each time-step") + xlab("Time") + 
    ylab("Number of Photons")

# get the estimated parameters
# estimated has a list of the estimated P_00, P_01, P_10, P_11, s0 an s1
  estimated = estimate_parameters(Y, 0.2, 0.3, 0.8, 600)

# run the functions many times to get different estimates with the same 
  #arguments
  sim_and_est = simulate_and_estimate(0.2, 0.3, 0.8, 50, 600, 200)
  #extract the estimates
  P_00_hat = sim_and_est$P_00_hat
  P_01_hat = sim_and_est$P_01_hat
  P_10_hat = sim_and_est$P_10_hat
  P_11_hat = sim_and_est$P_11_hat
  s0_hat = sim_and_est$s0_hat
  s1_hat=sim_and_est$s1_hat
  
  #plot in histograms the frequency of the estimates
  hist(P_00_hat,  main="Histogram of the estimated values of P_{00}",
       sub = "True P_{00} is equal to 0.3", xlab="Estimates for P_{00}", 
       border="yellow", labels = TRUE)
  
  hist(P_01_hat,  main="Histogram of the estimated values of P_{01}",
       sub = "True P_{01} is equal to 0.7", xlab="Estimates for P_{01}", 
       border="yellow")
  
  hist(P_10_hat,  main="Histogram of the estimated values of P_{10}", 
       sub = "True P_{10} is equal to 0.2", xlab="Estimates for P_{10}", 
       border="yellow")
  
  hist(P_11_hat,  main="Histogram of the estimated values of P_{11}", 
       sub = "True P_{11} is equal to 0.8", xlab="Estimates for P_{11}", 
       border="yellow")
  
  hist(s0_hat,  main="Histogram of the estimated values of s_{0}", 
       xlab="Estimates for s_{0}", 
       border="blue")
  
  hist(s1_hat,  main="Histogram of the estimated values of s_{1}", 
       xlab="Estimates for s_{1}", 
       border="blue")
  
  # calculate mse values between the predicted values and the actual ones
  mse_p00 = mse(0.3, P_00_hat)
  mse_p01 = mse(0.7, P_01_hat)
  mse_p10 = mse(0.2, P_10_hat)
  mse_p11 = mse(0.8, P_11_hat)

#### USE FUNCTIONS WITH NOISE
  
  Y = simulate_time_series_noise(X, 600, 200)
  
  p = estimate_parameters_noise(Y, 0.2, 0.3, 0.8, 600, 200)
  
  # plot the realization of photons in each step
  library(ggplot2)
  time = seq(0, 50)
  data_y_noise = cbind.data.frame(time, Y)
  
  ggplot(data_y_noise, aes(time, Y)) + geom_step() +
    ggtitle("Photon emissions in each time-step when using noise in the model") + xlab("Time") + 
    ylab("Number of Photons")
  
####Investigate the Baum-Welch algorithm
#Suppose that the realizations of a Markov chain are given, without having
# knowledge on the lambda and noise level. 
#Also, suppose that the transition probabilities are unknown as well.
# Try to estimate the transition probabilities in those cases.
  
  X_1 = simulate_MC(0.2, 0.3, 0.8, 250)
  library(ggplot2)
  time = seq(0, 250)
  data_x = cbind.data.frame(time, X_1)
  
  ggplot(data_x, aes(time, X_1)) + ylim(0,3) + geom_step() +
    ggtitle("Hidden States for each time-step") + xlab("Time") + 
    ylab("State")
  
  Y_1 = simulate_time_series_noise(X_1, 500, 350)
  library(ggplot2)
  time = seq(0, 250)
  data_y_noise = cbind.data.frame(time, Y_1)
  ggplot(data_y_noise, aes(time, Y_1)) + geom_step() +
    ggtitle("Photon emissions in each time-step, using noise in the model") + xlab("Time") + 
    ylab("Number of Photons")
  
  #Lambda and noise level unknown
  estimate_1 = estimate_parameters_noise(Y_1, 0.2, 0.3, 0.8, 500 ,350 )
  estimate_1$P_00_hat 
  estimate_1$P_01_hat 
  estimate_1$P_10_hat 
  estimate_1$P_11_hat 
  estimate_1$lambda_hat 
  
  #Lambda, noise level and transition probabilities unknown
  #Give as arguments the possible probabilities, lambda and noise based on the
  # figure which indicates the emitted photons and then estimate those values
  
  estimate_2 = estimate_parameters_noise(Y_1, 0.1, 0.5 ,0.9 , 700 ,250 )
  estimate_1$P_00_hat 
  estimate_1$P_01_hat 
  estimate_1$P_10_hat 
  estimate_1$P_11_hat 
  estimate_1$lambda_hat 
  
  estimate_3 = estimate_parameters_noise(Y_1, 0.3, 0.6 ,0.7 , 450 ,250 )
  estimate_3$P_00_hat
  estimate_3$P_01_hat 
  estimate_3$P_10_hat 
  estimate_3$P_11_hat 
  estimate_3$lambda_hat  
  