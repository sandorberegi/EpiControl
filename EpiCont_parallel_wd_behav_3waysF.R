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
  gen_time_var = c(2.1, 2.1),
  CFR = c(0.01, 0.5),
  mortality_mean = c(14.5, 10.0),
  mortality_var = c(1.1, 1.1)
)

behave_pars <- data.frame (
  epsilon = 0.35,
  cost_q = 20000,
  cost_diff = 1000,
  ind_impact = 0.0,
  cost_O = 0.25,
  eta2 = 0.9,
  softmax_coeff = 3000,
  cost_cum = 800000,
  cost_cum_shape = 12,
  emo_lam = 0.025,
  emo_R = 80000,
  emo_fatigue_coeff_N = 0.85,
  emo_fatigue_coeff_P = 0.7,
  emo_P_coeff = 0.25
)

r_trans_steep <- 1.5  # Growth rate
r_trans_len <- 7  # Number of days for the transition
t0 <- r_trans_len  / 2 # Midpoint of the transition

logistic_function <- function(t, R0, R1, r, t0) {
  K <- R1
  L <- R0
  (K-L) / (1 + exp(-r_trans_steep * (t - t0))) + L
}

I0 <- 10L #initial no. of infections
ndays <- 41L*7L#epidemic length
N <- 10e6 #Total population (if we account for susceptibles)

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
D_target <- 0
D_target_pen <- 500 #max death
alpha <- 1.3/C_target #~proportional gain (regulates error in cases) covid
#alpha = 3.25/C_target #~proportional gain (regulates error in cases) ebola
#beta <- 0.0 #~derivative gain (regulates error in R)
alpha_d <- 0*1.3/C_target*100
ovp <- 5.0 #overshoot penalty
dovp <- 0*10.0 #death overshoot penalty
gamma <- 0.95 #discounting factor

#Simulation parameters
n_ens <- 20L #MC assembly size for 4
sim_ens <- 20L #assembly size for full simulation

#Frequecy of policy review
rf <- 7L #days 14
R_est_wind <- 5L #rf-2 #window for R estimation
use_S <- 0L

#Prediction window
pred_days <- 12L #12L #14 #21 #12

# Original episim_data
column_names <- c("days", "sim_id", "I", "Lambda", "C", "Lambda_C", "S", "Deaths", "Re", "Rew", "Rest", "R0est", "policy", "R_coeff", "p", "n", "o", "pp", "Reserves", "LD_count", "SD_count", "uO_cost", "uP_cost", "uN_cost", "Cost_cum_O", "Cost_cum_P", "Cost_emo_N", "Cost_emo_P")

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
episim_data[1,] <- c(1, 1, I0, I0, Noise_pars['ur_mean']*I0, Noise_pars['ur_mean']*I0, N-I0, 0, Epi_pars[1,'R0'], Epi_pars[1,'R0'], 1, 1, 1, 1, 1, 0, 0, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0)

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
  episim_data_ens[[jj]] <- Epi_MPC_run_wd_behav_3waysF(episim_data_ens[[jj]], Epi_pars, Noise_pars, Action_space, behave_pars, pred_days = pred_days, start_day = 1, n_ens = n_ens, ndays = nrow(episim_data), R_est_wind = R_est_wind, pathogen = 1, susceptibles = 1, delay = 1, ur = 1, r_dir = 2, N = N)
}, cl=cl)

episim_data_ens <- results
stopCluster(cl)

for (jj in 1:sim_ens) {
  episim_data_ens[[jj]] <- head(episim_data_ens[[jj]], -pred_days)
}



#for (jj in 1:sim_ens) {
#  episim_data_ens[[jj]] <- Epi_MPC_run_wd_behav_3ways(episim_data_ens[[jj]], Epi_pars, Noise_pars, Action_space, behave_pars, pred_days = pred_days, start_day = 1, n_ens = n_ens, ndays = nrow(episim_data), R_est_wind = R_est_wind, pathogen = 1, susceptibles = 0, delay = 0, ur = 0, r_dir = 2, N = N)
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
  geom_line(data = subset(combined_data, sim_id != 1), aes(x = days, y = Deaths, color = as.factor(sim_id)), alpha = 0.3) +
  geom_line(data = subset(combined_data, sim_id == 1), aes(x = days, y = Deaths), color = "darkred", alpha = 1.0) +
  labs(x = "Days", y = "Deaths") +
  scale_color_manual(values = cls) +
  guides(color = FALSE)

ggplot() +
  geom_line(data = subset(combined_data, sim_id != 1), aes(x = days, y = C, color = as.factor(sim_id)), alpha = 0.3) +
  geom_line(data = subset(combined_data, sim_id == 1), aes(x = days, y = I), color = "darkred", alpha = 1.0) +
  geom_line(data = subset(combined_data, sim_id == 1), aes(x = days, y = C), color = "red", alpha = 1.0) +
  geom_line(data = subset(combined_data, sim_id == 1), aes(x = days, y = 200*Deaths), color = "blue", alpha = 1.0) +
  geom_hline(yintercept = C_target, linetype = "dashed", color = "purple", size=0.25) +
  geom_hline(yintercept = 200*D_target_pen, linetype = "dashed", color = "blue", size=0.25) +
  scale_y_continuous(
    name = "I and C",
    sec.axis = sec_axis(~./200, name = "Deaths")
  ) +
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
  geom_line(data = subset(combined_data, sim_id != 1), aes(x = days, y = Rest, color = as.factor(sim_id)), alpha = 0.3) +
  geom_line(data = subset(combined_data, sim_id == 1), aes(x = days, y = Re), color = "darkred", alpha = 1.0) +
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
combined_data


# Sample constant value
HI <- 1/Epi_pars[1,"R0"]

# Assuming episim_data_ens[[1]]["S"] is a dataframe or list
# Extract the vector
S_vector <- episim_data_ens[[1]]["S"][[1]]

# Compute the fraction S/N
fraction <- S_vector / N

# Find the index where the fraction is first smaller than the constant
index <- which(fraction < HI)[1]

# Output the result
index

# Plot with continuous lines and custom labels
ggplot(combined_data %>% filter(sim_id == 1)) +
  geom_line(data = subset(combined_data, sim_id != 1), aes(x = days, y = C, color = as.factor(sim_id)), alpha = 0.1) +
  geom_line(aes(x = days, y = C, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0) +
  geom_line(aes(x = days, y = 5000*n, color = "black" , group = 1), alpha = 1.0, size=0.25) +
  geom_line(aes(x = days, y = 5000*o, linetype = "dashed", color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0, size=0.25) +
  geom_hline(yintercept = C_target, linetype = "dashed", color = "blue", size=0.25) +
  geom_vline(xintercept = index, linetype = "dashed", color = "black", size=0.25) +
  scale_y_continuous(
    name = "Reported cases and true infections",
    sec.axis = sec_axis(~./5000, name = "p")
  ) +
  labs(x = "Days", y = "Reported cases and true infections", color = "Policy") +
  scale_color_manual(values = c("No intervention" = "chartreuse3", "Social distancing" = "darkorchid1", "Lockdown" = "red")) +
  guides(color = guide_legend(title = "Policy"))


# Plot with continuous lines and custom labels
ggplot(combined_data %>% filter(sim_id == 1)) +
  geom_point(data = subset(combined_data, sim_id != 1), aes(x = Re, y = p, color = as.factor(sim_id)), alpha = 0.1) +
  geom_point(aes(x = Re, y = p, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0) +
  labs(x = "Rt", y = "p", color = "Policy") +
  scale_color_manual(values = c("No intervention" = "chartreuse3", "Social distancing" = "darkorchid1", "Lockdown" = "red")) +
  guides(color = guide_legend(title = "Policy"))

# Plot with continuous lines and custom labels
ggplot(combined_data %>% filter(sim_id == 1)) +
  geom_point(data = subset(combined_data, sim_id != 1), aes(x = C, y = p, color = as.factor(sim_id)), alpha = 0.1) +
  geom_point(aes(x = C, y = p, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0) +
  labs(x = "Cases", y = "p", color = "Policy") +
  scale_color_manual(values = c("No intervention" = "chartreuse3", "Social distancing" = "darkorchid1", "Lockdown" = "red")) +
  guides(color = guide_legend(title = "Policy"))

ggplot(combined_data %>% filter(sim_id == 1)) +
  geom_point(data = subset(combined_data, sim_id != 1), aes(x = Lambda_C, y = p, color = as.factor(sim_id)), alpha = 0.1) +
  geom_point(aes(x = Lambda_C, y = p, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0) +
  labs(x = "Lambda_C", y = "p", color = "Policy") +
  scale_color_manual(values = c("No intervention" = "chartreuse3", "Social distancing" = "darkorchid1", "Lockdown" = "red")) +
  guides(color = guide_legend(title = "Policy"))

# Plot with continuous lines and custom labels
ggplot(combined_data %>% filter(sim_id == 1)) +
  geom_point(data = subset(combined_data, sim_id != 1), aes(x = C, y = Re, color = as.factor(sim_id)), alpha = 0.1) +
  geom_point(aes(x = C, y = Re, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0) +
  labs(x = "Cases", y = "R", color = "Policy") +
  scale_color_manual(values = c("No intervention" = "chartreuse3", "Social distancing" = "darkorchid1", "Lockdown" = "red")) +
  guides(color = guide_legend(title = "Policy"))

ggplot() +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = -uO_cost), color = "blue", alpha = 1.0) +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = -uP_cost), color = "darkgreen", alpha = 1.0) +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = -uN_cost), color = "red", alpha = 1.0) +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = -Cost_cum_O), color = "darkblue", alpha = 1.0) +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = -Cost_cum_P), color = "green", alpha = 1.0) +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = -Cost_emo_P), color = "darkgreen", alpha = 0.3) +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = -Cost_emo_N), color = "darkred", alpha = 1.0) +
  labs(x = "Days", y = "costs") +
  scale_color_manual(values = cls) +
  guides(color = FALSE)

ggplot() +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = -uO_cost), color = "blue", alpha = 1.0) +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = -uP_cost), color = "darkgreen", alpha = 1.0) +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = -uN_cost), color = "red", alpha = 1.0) +
  labs(x = "Days", y = "costs") +
  scale_color_manual(values = cls) +
  guides(color = FALSE)

ggplot() +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = o), color = "blue", alpha = 1.0) +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = p), color = "darkgreen", alpha = 1.0) +
  geom_step(data = subset(combined_data, sim_id == 1), aes(x = days, y = n), color = "red", alpha = 1.0) +
  labs(x = "Days", y = "costs") +
  scale_color_manual(values = cls) +
  guides(color = FALSE)

#filename <- paste0("behave_data__eps0.99_cost", behave_pars['cost_q'], ".csv")

# Write the dataframe to a CSV file
#write.csv(combined_data, filename, row.names = FALSE)

# Plot with continuous lines and custom labels
ggplot(combined_data %>% filter(sim_id == 1)) +
  geom_line(data = subset(combined_data, sim_id != 1), aes(x = days, y = C, color = as.factor(sim_id)), alpha = 0.1) +
  geom_line(aes(x = days, y = C, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0) +
  geom_line(aes(x = days, y = -Cost_emo_N/10, color = "black" , group = 1), alpha = 1.0, size=0.25) +
  scale_y_continuous(
    name = "Reported cases and true infections",
    sec.axis = sec_axis(~.*10, name = "Emotional cost")
  ) +
  labs(x = "Days", y = "Reported cases and true infections", color = "Policy") +
  scale_color_manual(values = c("No intervention" = "chartreuse3", "Social distancing" = "darkorchid1", "Lockdown" = "red")) +
  guides(color = guide_legend(title = "Policy"))
