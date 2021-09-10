X = simulate_fHMM(s0=c(0.2,0.4,0.6), P_00=c(0.2, 0.3, 0.8), P_11=c(0.7, 0.6, 0.5), t=50)

library(ggplot2)
time = seq(0, 50)
data_x = cbind.data.frame(time, t(X))

#plot of hidden states chains in separate geom_steps
ggplot() + ylim(0,3) + geom_step(data_x, mapping=aes(x=time, y=X[1,]), colour="red")+ geom_step(data_x, mapping=aes(x=time, y=X[2,]), colour="green")+
  geom_step(data_x, mapping=aes(x=time, y=X[3,]), colour="blue") + ggtitle("Hidden states for each molecule in each time-step", subtitle = "3 hidden chains in the FHMM") + xlab("Time-Step") + 
  ylab("Hidden States")

#plot of the hidden states chains
data_x = cbind.data.frame(time, colSums(X))
ggplot(data_x, mapping=aes(x=time, y=colSums(X))) + ylim(0, 5) + geom_step()+
  ggtitle("Number of molecules that are being activated in each time-step", subtitle = "3 hidden chains in the FHMM") +
  xlab("Time-Step") + ylab("Number of activated molecules")

#plot the emitted photons from fHMM, without noise
Y = simulate_time_series_fHMM(X, lambda=c(6000, 4000, 3000))
data_y = cbind.data.frame(time, Y)
ggplot(data_y, aes(time, Y)) + geom_line() +
  ggtitle("Photon emissions from all the molecules in each time-step", subtitle="Signal used: 6000, 4000 and 3000") + xlab("Time-Step") + 
  ylab("Number of emitted photons")

#plot the emitted photons from fHMM, with noise
Y_n = simulate_time_series_noise_fHMM(X, lambda=c(6000, 4000, 3000), lambda_n=2000)
data_y = cbind.data.frame(time, Y_n)
ggplot(data_y, aes(time, Y_n)) + geom_line() +
  ggtitle("Photon emissions from all the molecules in each time-step", subtitle ="Signal used: 6000, 4000 and 3000, Noise used: 2000") + xlab("Time-Step") + 
  ylab("Number of emitted photons")
