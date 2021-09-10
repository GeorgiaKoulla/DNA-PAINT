##Use the "estimate_frames_probs" function to simulate possible markov chains,
# given the average duration of the markov chain in the On and Off states.
# Based on that simulation, we estimate the emitted photons in each time-step.

# Getting possible simulation for the states of the markov chain when in slow, medium
# and fast switching scenarios

library(ggpubr)
#get the probabilities for each scenario
slow = estimate_frames_probs(0.2, 12.5)
medium = estimate_frames_probs(0.9, 16.7)
fast = estimate_frames_probs(3, 14.3)

#simulate markov chain states
s0 = 0.5 #set the initial probability that the chain starts at state 0, to be 0.5
X_slow = simulate_MC(s0, slow$p_00, slow$p_11, 3500)
X_medium = simulate_MC(s0, medium$p_00, medium$p_11, 3500)
X_fast = simulate_MC(s0, fast$p_00, fast$p_11, 3500)

#based on simulations estimate the photon emissions
Y_slow = simulate_time_series_noise(X_slow, lambda = 5000, 400)
Y_medium = simulate_time_series_noise(X_medium, lambda = 5000, 400)
Y_fast = simulate_time_series_noise(X_fast, lambda = 5000, 400)

#get the data ready
library(ggplot2)
time = seq(0, 3500)
data_slow = cbind.data.frame(time, X_slow)
data_medium = cbind.data.frame(time, X_medium)
data_fast = cbind.data.frame(time, X_fast)

data_yslow = cbind.data.frame(time, Y_slow)
data_ymedium = cbind.data.frame(time, Y_medium)
data_yfast = cbind.data.frame(time, Y_fast)

# SLOW PHOTO-SWITCHING SCENARIO
s1 = ggplot(data_slow, aes(time, X_slow)) + ylim(0,3) + geom_step() +
  ggtitle("States according to each time-step", 
          subtitle="Slow photo-switching scenario") + xlab("Time-Steps") + 
  ylab("State")
s2 = ggplot(data_yslow, aes(time, Y_slow)) + geom_step() +
  ggtitle("Photons that are emmited in each time-step",
          subtitle="Use signal 5000 and noise 400")+ xlab("Time-Steps") +
  ylab("Number of Photons")

fig1=ggarrange(s1, s2, ncol=2, nrow=1)
annotate_figure(fig1,top = text_grob("Slow photo-switching scenario", color = "green", face = "bold", size = 13))

# MEDIUM PHOTO-SWITCHING SCENARIO
m1 = ggplot(data_medium, aes(time, X_medium)) + ylim(0,3) + geom_step() +
  ggtitle("States according to each time-step", 
          subtitle="Medium photo-switching scenario") + xlab("Time-Steps") + 
  ylab("State")
m2 = ggplot(data_ymedium, aes(time, Y_medium)) + geom_step() +
  ggtitle("Photons that are emmited in each step",
          subtitle="Use signal 5000 and noise 400")+ xlab("Time-Steps") +
  ylab("Number of Photons")

fig2=ggarrange(m1, m2, ncol=2, nrow=1)
annotate_figure(fig2,top = text_grob("Medium photo-switching scenario", color = "blue", face = "bold", size = 13))

# FAST PHOTO-SWITCHING SCENARIO
f1 = ggplot(data_fast, aes(time, X_fast)) + ylim(0,3) + geom_step() +
  ggtitle("States according to each time-step", 
          subtitle="Fast photo-switching scenario") + xlab("Time-Steps") + 
  ylab("State")
f2 = ggplot(data_yfast, aes(time, Y_fast)) + geom_step() +
  ggtitle("Photons that are emmited in each step",
          subtitle="Use signal 5000 and noise 400")+ xlab("Time-Steps") +
  ylab("Number of Photons")

fig3=ggarrange(f1, f2, ncol=2, nrow=1)
annotate_figure(fig3,top = text_grob("Fast photo-switching scenario", color = "red", face = "bold", size = 13))


#####Estimate transition probabilities using the probabilities from each 
# photo-switching scenario.

#We will set the initial probability s0, that the chain starts at state 0 to be
# 0.5. So the probabilities that the chain starts from On and Off state are
# the same.

#SLOW PHOTO-SWITCHING CASE
s_case = smf_investigate_noise(0.2, 12.5, s0=0.5, t=1000, iter=100, lambda=5000,
                               noise_min=0, noise_max=2000)

s_mseP00 = s_case$mse_P_00_hat_level
s_mseP01 = s_case$mse_P_01_hat_level
s_mseP10 = s_case$mse_P_10_hat_level
s_mseP11 = s_case$mse_P_11_hat_level

s_averP00=s_case$aver_P_00_hat_level
s_averP01=s_case$aver_P_01_hat_level
s_averP10=s_case$aver_P_10_hat_level
s_averP11=s_case$aver_P_11_hat_level

s_realfrb = s_case$real_frb 
s_averfrb = s_case$aver_frb 

s_P00 = s_case$P_00
s_P11 = s_case$P_11
s_P01 = 1-s_P00
s_P10 = 1-s_P11

noise = seq(0, 2000, by=50)

#Average values for the estimated transition probabilities for each noise level

data=cbind.data.frame(noise, s_averP00)
s_a1 = ggplot(data=data, aes(noise, s_averP00))+ ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_00 for each noise level",
          subtitle = paste("Actual probability is", s_P00)) + 
  geom_hline(yintercept = s_P00, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_00")

data=cbind.data.frame(noise, s_averP01)
s_a2 = ggplot(data=data, aes(noise, s_averP01)) + ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_01 for each noise level", 
          subtitle = paste("Actual probability is", s_P01)) +
  geom_hline(yintercept = s_P01, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_01")

data=cbind.data.frame(noise, s_averP10)
s_a3 = ggplot(data=data, aes(noise, s_averP10)) + ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_10 for each noise level",
          subtitle=paste("Actual probability is", s_P10)) +
  geom_hline(yintercept = s_P10, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_10")

data=cbind.data.frame(noise, s_averP11)
s_a4 = ggplot(data=data, aes(noise, s_averP11)) + ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_11 for each noise level",
          subtitle =paste("Actual probability is", s_P11)) +
  geom_hline(yintercept = s_P11, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_11")

#MSE for the estimated transition probabilities for each noise level

data1=cbind.data.frame(noise, s_mseP00)
s_data1 = cbind.data.frame(data1, P_00 = "P_00")
names(s_data1)[2] = "MSE"
names(s_data1)[3] = "Transition_Probabilities"

data2=cbind.data.frame(noise, s_mseP01)
s_data2 = cbind.data.frame(data2, P_01 = "P_01")
names(s_data2)[2] = "MSE"
names(s_data2)[3] = "Transition_Probabilities"

data3=cbind.data.frame(noise, s_mseP10)
s_data3 = cbind.data.frame(data3, P_10 = "P_10")
names(s_data3)[2] = "MSE"
names(s_data3)[3] = "Transition_Probabilities"

data4=cbind.data.frame(noise, s_mseP11)
s_data4 = cbind.data.frame(data4, P_11 = "P_11")
names(s_data4)[2] = "MSE"
names(s_data4)[3] = "Transition_Probabilities"

df1 = rbind.data.frame(s_data1, s_data2, s_data3, s_data4) 

s_mse = ggplot(df1, aes(x = noise, y = MSE)) + ylim(0,0.75)+
  geom_line(aes(color = Transition_Probabilities, linetype = Transition_Probabilities)) + 
  scale_color_manual(values = c("red", "blue", "gold1","green"))+
  ggtitle("MSE for the estimates of the transition probabilities for each noise level")

# Average Frobenius norm of the error matrix for each noise level

data=cbind.data.frame(noise, s_averfrb)
s_f = ggplot(data=data, aes(noise, s_averfrb))+geom_line()+ylim(0, 1.75)+
  ggtitle("Average frobenius norm for each noise level", subtitle="Frobenius norm of the error matrix")+
  xlab("Noise Level") + ylab("Average frobenius norm")

s_fig=ggarrange(s_a1, s_a2, s_a3, s_a4, s_mse, s_f, ncol=2, nrow=3)
annotate_figure(s_fig,top = text_grob("Slow photo-switching scenario", color = "green",
                                         face = "bold", size = 13),
                fig.lab = "Red dashed line=actual probabilities", fig.lab.pos = "top.right", fig.lab.face="bold")


#MEDIUM PHOTO-SWITCHING CASE
m_case = smf_investigate_noise(0.9, 16.7, s0=0.5, t=1000, iter=100, lambda=5000,
                               noise_min=0, noise_max=2000)
m_mseP00 = m_case$mse_P_00_hat_level
m_mseP01 = m_case$mse_P_01_hat_level
m_mseP10 = m_case$mse_P_10_hat_level
m_mseP11 = m_case$mse_P_11_hat_level

m_averP00=m_case$aver_P_00_hat_level
m_averP01=m_case$aver_P_01_hat_level
m_averP10=m_case$aver_P_10_hat_level
m_averP11=m_case$aver_P_11_hat_level

m_realfrb = m_case$real_frb #1.376345 (l=700, noise_max=700)
m_averfrb = m_case$aver_frb 

m_P00 = m_case$P_00
m_P11 = m_case$P_11
m_P01 = 1-m_P00
m_P10 = 1-m_P11

noise = seq(0, 2000, by=50)

#Average values for the estimated transition probabilities for each noise level

data=cbind.data.frame(noise, m_averP00)
m_a1 = ggplot(data=data, aes(noise, m_averP00))+ ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_00 for each noise level",
          subtitle = paste("Actual probability is", m_P00)) + 
  geom_hline(yintercept = m_P00, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_00")

data=cbind.data.frame(noise, m_averP01)
m_a2 = ggplot(data=data, aes(noise, m_averP01)) + ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_01 for each noise level", 
          subtitle = paste("Actual probability is", m_P01)) +
  geom_hline(yintercept = m_P01, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_01")

data=cbind.data.frame(noise, m_averP10)
m_a3 = ggplot(data=data, aes(noise, m_averP10)) + ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_10 for each noise level",
          subtitle=paste("Actual probability is", m_P10)) +
  geom_hline(yintercept = m_P10, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_10")

data=cbind.data.frame(noise, m_averP11)
m_a4 = ggplot(data=data, aes(noise, m_averP11)) + ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_11 for each noise level",
          subtitle =paste("Actual probability is", m_P11)) +
  geom_hline(yintercept = m_P11, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_11")

#MSE for the estimated transition probabilities for each noise level

data1=cbind.data.frame(noise, m_mseP00)
m_data1 = cbind.data.frame(data1, P_00 = "P_00")
names(m_data1)[2] = "MSE"
names(m_data1)[3] = "Transition_Probabilities"

data2=cbind.data.frame(noise, m_mseP01)
m_data2 = cbind.data.frame(data2, P_01 = "P_01")
names(m_data2)[2] = "MSE"
names(m_data2)[3] = "Transition_Probabilities"

data3=cbind.data.frame(noise, m_mseP10)
m_data3 = cbind.data.frame(data3, P_10 = "P_10")
names(m_data3)[2] = "MSE"
names(m_data3)[3] = "Transition_Probabilities"

data4=cbind.data.frame(noise, m_mseP11)
m_data4 = cbind.data.frame(data4, P_11 = "P_11")
names(m_data4)[2] = "MSE"
names(m_data4)[3] = "Transition_Probabilities"

df2 = rbind.data.frame(m_data1, m_data2, m_data3, m_data4) 

m_mse = ggplot(df2, aes(x = noise, y = MSE)) + ylim(0,0.75)+
  geom_line(aes(color = Transition_Probabilities, linetype = Transition_Probabilities)) + 
  scale_color_manual(values = c("red", "blue", "gold1","green"))+
  ggtitle("MSE for the estimates of the transition probabilities for each noise level")

# Average Frobenius norm of the error matrix for each noise level

data=cbind.data.frame(noise, m_averfrb)
m_f = ggplot(data=data, aes(noise, m_averfrb))+geom_line()+ylim(0, 1.75)+
  ggtitle("Average frobenius norm for each noise level", subtitle="Frobenius norm of the error matrix")+
  xlab("Noise Level") + ylab("Average frobenius norm")

m_fig=ggarrange(m_a1, m_a2, m_a3, m_a4, m_mse, m_f, ncol=2, nrow=3)
annotate_figure(m_fig, top = text_grob("Medium photo-swithcing scenario", color = "blue",
                                      face = "bold", size = 13),
                fig.lab = "Red dashed line=actual probabilities", fig.lab.pos = "top.right", fig.lab.face="bold")


#FAST PHOTO-SWITCHING CASE
f_case = smf_investigate_noise(3, 14.3, s0=0.5, t=1000, iter=100, lambda=5000,
                               noise_min=0, noise_max=2000)

f_mseP00 = f_case$mse_P_00_hat_level
f_mseP01 = f_case$mse_P_01_hat_level
f_mseP10 = f_case$mse_P_10_hat_level
f_mseP11 = f_case$mse_P_11_hat_level

f_averP00=f_case$aver_P_00_hat_level
f_averP01=f_case$aver_P_01_hat_level
f_averP10=f_case$aver_P_10_hat_level
f_averP11=f_case$aver_P_11_hat_level

f_realfrb = f_case$real_frb 
f_averfrb = f_case$aver_frb 

f_P00 = f_case$P_00
f_P11 = f_case$P_11
f_P01 = 1-f_P00
f_P10 = 1-f_P11

noise = seq(0, 2000, by=50)

#Average values for the estimate transition probabilities for each noise level

data=cbind.data.frame(noise, f_averP00)
f_a1 = ggplot(data=data, aes(noise, f_averP00))+ ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_00 for each noise level",
          subtitle = paste("Actual probability is", f_P00)) + 
  geom_hline(yintercept = f_P00, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_00")

data=cbind.data.frame(noise, f_averP01)
f_a2 = ggplot(data=data, aes(noise, f_averP01)) + ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_01 for each noise level", 
          subtitle = paste("Actual probability is", f_P01)) +
  geom_hline(yintercept = f_P01, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_01")

data=cbind.data.frame(noise, f_averP10)
f_a3 = ggplot(data=data, aes(noise, f_averP10)) + ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_10 for each noise level",
          subtitle=paste("Actual probability is", f_P10)) +
  geom_hline(yintercept = f_P10, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_10")

data=cbind.data.frame(noise, f_averP11)
f_a4 = ggplot(data=data, aes(noise, f_averP11)) + ylim(0,2) + geom_line()+
  ggtitle("Average values for the estimates of P_11 for each noise level",
          subtitle =paste("Actual probability is", f_P11)) +
  geom_hline(yintercept = f_P11, linetype="dashed", color="red")+
  xlab("Noise Level") + ylab("Average values for the estimates of P_11")

#MSE for the estimated transition probabilities for each noise level

data1=cbind.data.frame(noise, f_mseP00)
f_data1 = cbind.data.frame(data1, P_00 = "P_00")
names(f_data1)[2] = "MSE"
names(f_data1)[3] = "Transition_Probabilities"

data2=cbind.data.frame(noise, f_mseP01)
f_data2 = cbind.data.frame(data2, P_01 = "P_01")
names(f_data2)[2] = "MSE"
names(f_data2)[3] = "Transition_Probabilities"

data3=cbind.data.frame(noise, f_mseP10)
f_data3 = cbind.data.frame(data3, P_10 = "P_10")
names(f_data3)[2] = "MSE"
names(f_data3)[3] = "Transition_Probabilities"

data4=cbind.data.frame(noise, f_mseP11)
f_data4 = cbind.data.frame(data4, P_11 = "P_11")
names(f_data4)[2] = "MSE"
names(f_data4)[3] = "Transition_Probabilities"

df3 = rbind.data.frame(f_data1, f_data2, f_data3, f_data4) 

f_mse = ggplot(df3, aes(x = noise, y = MSE)) + ylim(0,0.75)+
  geom_line(aes(color = Transition_Probabilities, linetype = Transition_Probabilities)) + 
  scale_color_manual(values = c("red", "blue", "gold1","green"))+
  ggtitle("MSE for the estimates of the transition probabilities for each noise level")

# Average Frobenius norm of the error matrix for each noise level

data=cbind.data.frame(noise, f_averfrb)
f_f = ggplot(data=data, aes(noise, f_averfrb))+geom_line()+ylim(0, 1.75)+
  ggtitle("Average frobenius norm for each noise level", subtitle="Frobenius norm of the error matrix")+
  xlab("Noise Level") + ylab("Average frobenius norm")

f_fig=ggarrange(f_a1, f_a2, f_a3, f_a4, f_mse, f_f, ncol=2, nrow=3)
annotate_figure(f_fig,top = text_grob("Fast photo-switching scenario", color = "red",
                                      face = "bold", size = 13),
                fig.lab = "Red dashed line=actual probabilities", fig.lab.pos = "top.right", fig.lab.face="bold")
