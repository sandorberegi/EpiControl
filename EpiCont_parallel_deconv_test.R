## Declaring epidemic & noise parameters

library(VGAM)
#library(future.apply)
#library(foreach)
library(parallel)
library(pbapply)

cores=detectCores()-1
cl <- makeCluster(cores)

Epi_pars = data.frame (
  Pathogen = c("COVID-19", "Ebola"),
  R0 = c(3.5, 2.5),
  gen_time = c(6.5, 15.0),
  gen_time_var = c(2.1, 2.1)
)

I0 <- 10L #initial no. of infections
ndays <- 51L*7L#epidemic length
N <- 1e7 #Total population (if we account for susceptibles)

Noise_pars <- data.frame (
  repd_mean = 10.5, #Reporting delay mean
  del_disp = 5.0, #Reporting delay variance
  ur_mean = 0.3, #Under reporting, mean/variance
  ur_beta_a = 50.0
)

# Setting-up the control (using non-pharmaceutical interventions)
Action_space <- data.frame (
  NPI = c("No restrictions", "Social distancing", "Lockdown"),
  R_coeff = c(1.0, 0.5, 0.2), #R0_act = R0 * ctrl_states
  R_beta_a = c(0.0, 5.0, 5.0), #R0_act uncertainty
  cost_of_NPI = c(0.0, 0.01, 0.15)
)

#Sim and control options
cost_sel <- 1L # : bilinear+oveshoot, 2: flat+quadratic for overshoot, 3: ReLu + Overshoot
use_inc <- 1L #1: control for incidence, 0: control for infectiousness
delay_calc_v <- 0L #1: as in ref, 0: from incidence
under_rep_calc <- 3L #1: as in ref, 0: separately, from incidence, 2: same, but using a beta distribution for rho
distr_sel <- 1L #0: Deterministic, 1: Poisson, 2: Binomial
delay_on <- 1L #1: sim with time-delay, 0: no-delay (if 0, set delay_calc_v = 0)
under_rep_on <- 1L #0: no under reporting, 1: calculate with under-reporting

C_target <- 5000 #target cases
C_target_pen <- C_target*1.5 #overshoot penalty threshold
R_target <- 1.0
alpha <- 1.3/C_target #~proportional gain (regulates error in cases) covid
#alpha = 3.25/C_target #~proportional gain (regulates error in cases) ebola
#beta <- 0.0 #~derivative gain (regulates error in R)
ovp <- 5.0 #overshoot penalty
gamma <- 0.95 #discounting factor

#Simulation parameters
n_ens <- 100L #MC assembly size for 4
sim_ens <- 100L #assembly size for full simulation

#Frequecy of policy review
rf <- 7L #days 14
R_est_wind <- 5L #rf-2 #window for R estimation
use_S <- 0L

#Prediction window
pred_days <- 12L #14 #21 #12

# Original episim_data
column_names <- c("days", "sim_id", "I", "Lambda", "C", "Lambda_C", "S", "Re", "Rew", "Rest", "R0est", "policy", "R_coeff")

# Create an empty data frame with specified column names
empty_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(empty_df) <- column_names

# Print the structure of the empty data frame
# Create a zero matrix with 10 rows
zero_matrix <- matrix(0, nrow = ndays, ncol = length(column_names))
colnames(zero_matrix) <- column_names

# Combine empty data frame and zero matrix using rbind
episim_data <- rbind(empty_df, zero_matrix)

#initialisation

episim_data['policy'] <- rep(1, ndays)
episim_data['sim_id'] <- rep(1, ndays)
episim_data['days'] <- 1:ndays
episim_data[1,] <- c(1, 1, I0, I0, Noise_pars['ur_mean']*I0, Noise_pars['ur_mean']*I0, N-I0, Epi_pars[1,'R0'], Epi_pars[1,'R0'], 1, 1, 1)

episim_data_ens <- replicate(sim_ens, episim_data, simplify = FALSE)

pb <- txtProgressBar(min = 0, max = sim_ens, initial = 0, style = 3, width = 50, char = "=")
for (ii in 1:sim_ens) {
  episim_data_ens[[ii]]$sim_id <- rep(ii, ndays)
}

# Parallel setup
#num_cores <- detectCores() - 1  # Use one less than the total number of cores
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)

#plan(multisession) # Use multisession plan for parallel execution

# Initialize the progress bar
pb <- txtProgressBar(min = 0, max = sim_ens, initial = 0, style = 3, width = 50, char = "=")

# Ensure each simulation has a unique sim_id
for (ii in 1:sim_ens) {
  episim_data_ens[[ii]]$sim_id <- rep(ii, ndays)
}

clusterExport(cl, ls())

clusterEvalQ(cl, episim_data_ens)
clusterEvalQ(cl, {
  library(VGAM)
  library(EpiControl)
})

results <- pblapply(1:sim_ens, function(jj) {
  episim_data_ens[[jj]] <- Epi_MPC_run(episim_data_ens[[jj]], Epi_pars, Noise_pars, Action_space, pred_days = pred_days, start_day = 1, n_ens = n_ens, ndays = nrow(episim_data), R_est_wind = R_est_wind, pathogen = 1, susceptibles = 0, delay = 1, ur = 0, N = N)
}, cl=cl)

episim_data_ens <- results
stopCluster(cl)

trial_data <- episim_data_ens[[1]]

repd_mean <- Noise_pars[1,'repd_mean'] #Reporting delay mean
del_disp <- Noise_pars[1,'del_disp'] #Reporting delay variance

Ydel <- dgamma(1:ndays, del_disp, del_disp/repd_mean)
Ydel <- (Ydel/sum(Ydel))[1:31]

plot(Ydel)

library(fastbeta)

set.seed(2L)
#n <- 200L
#d <- 50L
#p <- 0.1
#prob <- plogis(rlogis(d + n, location = qlogis(p), scale = 0.1))


#delay <- diff(pgamma(0L:(d + 1L), 12, 0.4))

#h <- function (x, a = 1, b = 1, c = 0) a * exp(-b * (x - c)^2)
#ans <- floor(h(seq(-60, 60, length.out = d + n), a = 1000, b = 0.001))

#x0 <- rbinom(d + n, ans, prob)
#x <- tabulate(rep.int(1L:(d + n), x0) +
#                sample(0L:d, size = sum(x0), replace = TRUE, prob = delay),
#              d + n)[-(1L:d)]

x <- trial_data[['C']]
y <- trial_data[['I']]

str(D0 <- deconvolve(x, 1.0, delay = Ydel, iter.max = 64L, complete = FALSE))
#str(D1 <- deconvolve(x, 1.0, delay = Ydel, complete =  TRUE))

plot(D0$value)
plot(x)

trial_data['Deconv_C'] <- D0$value

# Plot x first
plot(x, type = "l", col = "blue", ylim = range(c(x, D0$value)), ylab = "Values", xlab = "Index", main = "Overlay of x and D0$value")

# Add D0$value to the same plot
lines((D0$value)[31:(387-15)], col = "red")
lines(y, col = "black")

# Optionally, add a legend to differentiate between the two
legend("topright", legend = c("x", "D0$value"), col = c("blue", "red"), lty = 1)









#for (jj in 1:sim_ens) {
#  episim_data_ens[[jj]] <- Epi_MPC_run(episim_data_ens[[jj]], Epi_pars, Noise_pars, Action_space, pred_days = pred_days, n_ens = n_ens, ndays = nrow(episim_data), R_est_wind = R_est_wind, pathogen = 1, susceptibles = 0, delay = 0, ur = 0, N = N)
#  setTxtProgressBar(pb,jj)
#}
#close(pb)

combined_data <- do.call(rbind, episim_data_ens)

# Print the combined data frame
print(combined_data)

library(ggplot2)

# Plot the ensemble data for columns 'days' and 'I'
ggplot(combined_data, aes(x = days, y = C)) +
  geom_line() +
  labs(title = "Ensemble Data for 'days' and 'C'",
       x = "Days",
       y = "Cases")

quantile_vals <- c(0.025, 0.5, 0.975)

library(dplyr)
summary_data <- combined_data %>%
  group_by(days) %>%
  summarise(q5 = quantile(C, probs = 0.05),
            mean = mean(C),
            q95 = quantile(C, probs = 0.95))

summary_data

ggplot() +
  geom_line(data = combined_data, aes(x = days, y = C), color = "grey", alpha = 0.3) +
  geom_line(data = summary_data, aes(x = days, y = mean), color = "red", size = 1) +
  geom_ribbon(data = summary_data, aes(x = days, ymin = q5, ymax = q95), alpha = 0.2, fill = "red") +
  labs(x = "Time", y = "Ensemble average")

cls <- rep("grey", sim_ens)
cls[1] <- "red"

ggplot() +
  geom_line(data = subset(combined_data, sim_id != 1), aes(x = days, y = C, color = as.factor(sim_id)), alpha = 0.3) +
  geom_line(data = subset(combined_data, sim_id == 1), aes(x = days, y = I), color = "darkred", alpha = 1.0) +
  geom_line(data = subset(combined_data, sim_id == 1), aes(x = days, y = C), color = "red", alpha = 1.0) +
  labs(x = "Days", y = "I") +
  scale_color_manual(values = cls) +
  guides(color = FALSE)

ggplot() +
  geom_line(data = subset(combined_data, sim_id != 1), aes(x = days, y = Rest, color = as.factor(sim_id)), alpha = 0.3) +
  geom_line(data = subset(combined_data, sim_id == 1), aes(x = days, y = Rew), color = "darkred", alpha = 1.0) +
  geom_line(data = subset(combined_data, sim_id == 1), aes(x = days, y = Rest), color = "red", alpha = 1.0) +
  labs(x = "Days", y = "R") +
  scale_color_manual(values = cls) +
  guides(color = FALSE)

ggplot() +
  geom_line(data = subset(combined_data, sim_id != 1), aes(x = days, y = R0est, color = as.factor(sim_id)), alpha = 0.3) +
  geom_line(data = subset(combined_data, sim_id == 1), aes(x = days, y = R0est), color = "red", alpha = 1.0) +
  labs(x = "Days", y = "R0 estimate") +
  scale_color_manual(values = cls) +
  guides(color = FALSE)

# Create a new grouping variable for discontinuities
combined_data <- combined_data %>%
  group_by(sim_id, policy) %>%
  arrange(days) %>%
  mutate(group = cumsum(c(1, diff(days) != 1))) %>%
  ungroup()

combined_data <- combined_data %>%
  arrange(sim_id, days, policy)

policy_labels <- c("1" = "No intervention", "2" = "Social distancing", "3" = "Lockdown")

# Plot with continuous lines and custom labels
ggplot(combined_data %>% filter(sim_id == 1)) +
  geom_line(data = subset(combined_data, sim_id != 1), aes(x = days, y = C, color = as.factor(sim_id)), alpha = 0.1) +
  geom_line(aes(x = days, y = C, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0) +
  geom_line(aes(x = days, y = I, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0, size=0.25) +
  geom_hline(yintercept = C_target, linetype = "dashed", color = "blue", size=0.25) +
  labs(x = "Days", y = "Reported cases and true infections", color = "Policy") +
  scale_color_manual(values = c("No intervention" = "chartreuse3", "Social distancing" = "darkorchid1", "Lockdown" = "red")) +
  guides(color = guide_legend(title = "Policy"))
