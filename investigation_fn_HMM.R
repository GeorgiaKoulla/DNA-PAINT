### INVESTIGATIONS
library(ggplot2)

#Setting s0, P_00 and P_11 to be small
X1 = simulate_MC(0.3, 0.2, 0.3, 100)
time = seq(0, 100)
data_x = cbind.data.frame(time, X1)
ggplot(data_x, aes(time, X1)) + ylim(0,3) + geom_step() +
  ggtitle("States for each t(time)") + xlab("Time") + 
  ylab("State")

#Setting s0, P_00 to be small and P_11 to be high
X2 = simulate_MC(0.3, 0.3, 0.7, 100 )
time = seq(0, 100)
data_x = cbind.data.frame(time, X2)
ggplot(data_x, aes(time, X2)) + ylim(0,3) + geom_step() +
  ggtitle("States for each t(time)") + xlab("Time") + 
  ylab("State")

#Setting s0, P_11 to be small and P_00 to be high
X3 = simulate_MC(0.3, 0.7, 0.3, 100 )
time = seq(0, 100)
data_x = cbind.data.frame(time, X3)
ggplot(data_x, aes(time, X3)) + ylim(0,3) + geom_step() +
  ggtitle("States for each t(time)") + xlab("Time") + 
  ylab("State")

#Setting s0 to be small and P_00 and P_11 to be high
X4 = simulate_MC(0.3, 0.8, 0.8, 100 )
time = seq(0, 100)
data_x = cbind.data.frame(time, X4)
ggplot(data_x, aes(time, X4)) + ylim(0,3) + geom_step() +
  ggtitle("States for each t(time)") + xlab("Time") + 
  ylab("State")

#Setting s0, P_00 and P_11 to be high
X5 = simulate_MC(0.7, 0.8, 0.6, 100 )
time = seq(0, 100)
data_x = cbind.data.frame(time, X5)
ggplot(data_x, aes(time, X5)) + ylim(0,3) + geom_step() +
  ggtitle("States for each t(time)") + xlab("Time") + 
  ylab("State")

#Setting s0, P_00 to be high and P_11 to be small
X6 = simulate_MC(0.7, 0.7, 0.4, 100 )
time = seq(0, 100)
data_x = cbind.data.frame(time, X6)
ggplot(data_x, aes(time, X6)) + ylim(0,3) + geom_step() +
  ggtitle("States for each t(time)") + xlab("Time") + 
  ylab("State")

#Setting s0, P_11 to be high and P_00 to be small
X7 = simulate_MC(0.7, 0.3, 0.8, 100 )
time = seq(0, 100)
data_x = cbind.data.frame(time, X7)
ggplot(data_x, aes(time, X7)) + ylim(0,3) + geom_step() +
  ggtitle("States for each t(time)") + xlab("Time") + 
  ylab("State")

#Setting s0 to be high and P_00 and P_11 to be small
X8 = simulate_MC(0.8, 0.2, 0.3, 100 )
time = seq(0, 100)
data_x = cbind.data.frame(time, X8)
ggplot(data_x, aes(time, X8)) + ylim(0,3) + geom_step() +
  ggtitle("States for each t(time)") + xlab("Time") + 
  ylab("State")

# INVESTIGATION OF THE REALIZATIONS(Y)-THE EMMITTED PHOTONS IN EACH STEP
 
# generate an X - states of the Markov chain in each time-step, according to 
# the probabilities s0, P_00, P_11 and then get the realizations of the photons

X_y = simulate_MC(0.3, 0.4, 0.7, 100)
time = seq(0, 100)
data_x = cbind.data.frame(time, X_y)
ggplot(data_x, aes(time, X_y)) + ylim(0,3) + geom_step() +
  ggtitle("Hidden States for each time-step") + xlab("Time") + 
  ylab("State")

Y_s = simulate_time_series(X_y, 600)
data_y = cbind.data.frame(time, Y_s)
ggplot(data_y, aes(time, Y_s)) + geom_step() +
  ggtitle("Photon Emissions in each time-step") + xlab("Time") + 
  ylab("Number of Photons")

Y_s_n = simulate_time_series_noise(X_y, 600, 300)
data_y = cbind.data.frame(time, Y_s_n)
ggplot(data_y, aes(time, Y_s_n)) + geom_step() +
  ggtitle("Photon Emissions in each time-step, using noise in the model") + xlab("Time") + 
  ylab("Number of Photons")

#INVESTIGATION OF THE MSE OF TRANSITION PROBABILITIES AND COMPARISONS OF
# TRANSITION MATRICES
estimations = simulate_and_estimate(0.3, 0.4, 0.7, 100, 600, 150)

P_00_hat = estimations$P_00_hat
P_01_hat = estimations$P_01_hat
P_10_hat = estimations$P_10_hat
P_11_hat = estimations$P_11_hat
s0_hat = estimations$s0_hat
s1_hat=estimations$s1_hat

mse_p00 = mse(0.4, P_00_hat) 
mse_p01 = mse(0.6, P_01_hat) 
mse_p10 = mse(0.3, P_10_hat) 
mse_p11 = mse(0.7, P_11_hat) 

estimations_n = simulate_and_estimate_noise(0.3, 0.4, 0.7, 100, 600, 200, 150)

P_00_hat_n = estimations_n$P_00_hat
P_01_hat_n = estimations_n$P_01_hat
P_10_hat_n = estimations_n$P_10_hat
P_11_hat_n = estimations_n$P_11_hat
s0_hat_n = estimations_n$s0_hat
s1_hat_n=estimations_n$s1_hat

mse_p00_n = mse(0.4, P_00_hat_n) 
mse_p01_n = mse(0.6, P_01_hat_n) 
mse_p10_n = mse(0.3, P_10_hat_n) 
mse_p11_n = mse(0.7, P_11_hat_n) 


#Comparison of the estimated transition matrices using the Frobenius method.

#frobenius norm of the transition matrix with the true values that we set
frob_norm(0.4, 0.7) #1.0488

#Frobenius norm of the transition matrix with the estimated probabilities
frob_norm(P_00_hat, P_11_hat) 

frob_norm(P_00_hat_n, P_11_hat_n) 

###INVESTIGATION BASED ON THE LEVEL OF THE NOISE

#Inspect the affect of the noise level when estimating the transition
# probabilities of the Markov chain

X_n = simulate_MC(0.4, 0.4, 0.7, 200)
library(ggplot2)
time = seq(0, 200)
data_x = cbind.data.frame(time, X_n)
ggplot(data_x, aes(time, X_n)) + ylim(0,3) + geom_step() +
  ggtitle("Hidden States for each time-step") + xlab("Time") + 
  ylab("State")

Y_n1 = simulate_time_series_noise(X_n, lambda = 600, 250)
data_y1 = cbind.data.frame(time, Y_n1)
Y_n2 = simulate_time_series_noise(X_n, lambda = 600, 300)
data_y2 = cbind.data.frame(time, Y_n2)
Y_n3 = simulate_time_series_noise(X_n, lambda = 600, 400)
data_y3 = cbind.data.frame(time, Y_n3)
Y_n4 = simulate_time_series_noise(X_n, lambda = 600, 500)
data_y4 = cbind.data.frame(time, Y_n4)

p1 = ggplot(data_y1, aes(time, Y_n1)) + geom_step() +
  ggtitle("Photon Emissions in each time-step, using noise level=250") + xlab("Time") + 
  ylab("Number of Photons")
p2 = ggplot(data_y2, aes(time, Y_n2)) + geom_step() +
  ggtitle("Photon Emissions in each time-step, using noise level=300") + xlab("Time") + 
  ylab("Number of Photons")
p3 = ggplot(data_y3, aes(time, Y_n3)) + geom_step() +
  ggtitle("Photon Emissions in each time-step, using noise level=400") + xlab("Time") + 
  ylab("Number of Photons")
p4 = ggplot(data_y4, aes(time, Y_n4)) + geom_step() +
  ggtitle("Photon Emissions in each time-step, using noise level=500") + xlab("Time") + 
  ylab("Number of Photons")

library(patchwork)
(p1 / p2)
(p3 /p4)
  
#set the model with s0=0.3, P_00=0.6, P_11=0.8, lambda=600, noise=0,50,...,600
noise_inv = investigate_noise_level(s0=0.3, P_00=0.6, P_11=0.8, lambda=600, 
                                    t=1000, iter=100, noise_min=0, noise_max=600)

mseP00 = noise_inv$mse_P_00_hat_level
mseP01 = noise_inv$mse_P_01_hat_level
mseP10 = noise_inv$mse_P_10_hat_level
mseP11 = noise_inv$mse_P_11_hat_level


averP00=noise_inv$aver_P_00_hat_level
averP01=noise_inv$aver_P_01_hat_level
averP10=noise_inv$aver_P_10_hat_level
averP11=noise_inv$aver_P_11_hat_level


realfrb = noise_inv$real_frb 
averfrb = noise_inv$aver_frb 

####NOW, INVESTIGATE THE ABOVE WITH PLOTS
noise = seq(0,600,by=50)

#Average values of the estimated transition probabilities for each noise level

data=cbind.data.frame(noise, averP00)
a1 = ggplot(data=data, aes(noise, averP00))+ ylim(0,2) + geom_line()+
  ggtitle("Average value for the estimates of P_00 for each noise level",
          subtitle ="Actual probability is 0.6") + 
  geom_hline(yintercept = 0.6, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_00")

data=cbind.data.frame(noise, averP01)
a2 = ggplot(data=data, aes(noise, averP01)) + ylim(0,2) + geom_line()+
  ggtitle("Average value for the estimates of P_01 for each noise level", 
          subtitle = "Actual probability is 0.4") +
  geom_hline(yintercept = 0.4, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_01")

data=cbind.data.frame(noise, averP10)
a3 = ggplot(data=data, aes(noise, averP10)) + ylim(0,2) + geom_line()+
  ggtitle("Average  value for the estimates of P_10 for each noise level",
          subtitle="Actual probability is 0.2") +
  geom_hline(yintercept = 0.2, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_10")

data=cbind.data.frame(noise, averP11)
a4 = ggplot(data=data, aes(noise, averP11)) + ylim(0,2) + geom_line()+
  ggtitle("Average value for the estimates of P_11 for each noise level",
          subtitle ="Actual probability is 0.8") +
  geom_hline(yintercept = 0.8, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_11")

aver_fig=ggarrange(a1, a2, a3, a4, ncol=2, nrow=2)
annotate_figure(aver_fig,top = text_grob("Average values for the estimates of the transition probabilities", color = "blue",
                                         face = "bold", size = 13),
                fig.lab = "Red dashed line=actual probabilities", fig.lab.pos = "top.right", fig.lab.face="bold")


#MSE for the estimated transition probabilities for each noise level

data=cbind.data.frame(noise, mseP00)
e1 = ggplot(data=data, aes(noise, mseP00))+geom_line()+ ylim(0,0.005)+
  ggtitle("MSE for the estimates of P_00 for each noise level")+
  xlab("Noise Level") + ylab("MSE for the estimates of P_00")

data=cbind.data.frame(noise, mseP01)
e2 = ggplot(data=data, aes(noise, mseP01))+geom_line()+ylim(0,0.005)+
  ggtitle("MSE for the estimates of P_01 for each noise level")+
  xlab("Noise Level") + ylab("MSE for the estimates of P_01")

data=cbind.data.frame(noise, mseP10)
e3 = ggplot(data=data, aes(noise, mseP10))+geom_line()+ylim(0,0.005)+
  ggtitle("MSE for the estimates of P_10 for each noise level")+
  xlab("Noise Level") + ylab("MSE for the estimates of P_10")

data=cbind.data.frame(noise, mseP11)
e4 = ggplot(data=data, aes(noise, mseP11))+geom_line()+ylim(0,0.005)+
  ggtitle("MSE for the estimates of P_11 for each noise level")+
  xlab("Noise Level") + ylab("MSE for the estimates of P_11")

mse_fig=ggarrange(e1, e2, e3, e4, ncol=2, nrow=2)
annotate_figure(mse_fig,top = text_grob("MSE for the estimates of the transition probabilities", color = "blue",
                                        face = "bold", size = 13))


#Average Frobenius norm for each noise level, from the difference of the 
# estimated transition matrix and the real transition matrix

data=cbind.data.frame(noise, averfrb)
f = ggplot(data=data, aes(noise, averfrb))+geom_line()+ylim(1.25, 1.75)+
  ggtitle("Average Frobenius norm for each noise level", subtitle="Frobenius norm of the error matrix: difference between estimated and real transition matrix")+
  xlab("Noise Level") + ylab("Average Frobenius norm")
