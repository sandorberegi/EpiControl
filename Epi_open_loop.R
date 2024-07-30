## Declaring epidemic & noise parameters

library(VGAM)

Epi_pars = data.frame (
  Pathogen = c("COVID-19", "Ebola"),
  R0 = c(3.5, 2.5),
  gen_time = c(6.5, 15.0),
  gen_time_var = c(2.1, 2.1)
)

I0 <- 10L #initial no. of infections
ndays <- 41L*7L#epidemic length
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
beta <- 0.0 #~derivative gain (regulates error in R)
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

for (ii in 1:sim_ens) {
  episim_data_ens[[ii]]$sim_id <- rep(ii, ndays)
}

#Epi_run <- function(episimdata, epi_par, noise_par) { # create a function with the name my_function
#
#}

# Epi_run <- function(episimdata, epi_par, noise_par, ndays = nrow(episimdata), pathogen = 1, susceptibles = 1, delay = 0, ur = 0) {
#
#   R0 <- epi_par[pathogen,"R0"]
#   gen_time <- epi_par[pathogen,"gen_time"]
#   gen_time_var <- epi_par[pathogen,"gen_time_var"]
#
#   Ygen <- dgamma(1:ndays, gen_time/gen_time_var, 1/gen_time_var)
#   Ygen <- Ygen/sum(Ygen)
#
#   if (delay == 1) {
#     repd_mean <- noise_par[1,'repd_mean'] #Reporting delay mean
#     del_disp <- noise_par[1,'del_disp'] #Reporting delay variance
#
#     Ydel <- dgamma(1:ndays, del_disp, del_disp/repd_mean)
#     Ydel <- Ydel/sum(Ydel)
#   }
#
#   if (ur == 1) {
#     ur_mean <- noise_par[1,'ur_mean'] #Under reporting, mean/variance
#     ur_beta_a <- noise_par[1,'ur_beta_a']
#     ur_beta_b <- (1-ur_mean)/ur_mean * ur_beta_a
#   }
#
#   for (ii in 2:ndays) {
#
#     if (susceptibles == 1) {
#       Ract <- R0 * episimdata[ii-1,'S']/N
#     } else {
#       Ract <- R0
#     }
#
#     episimdata[ii, 'Re'] <- Ract
#     episimdata[ii, 'Rew'] <- sum(episimdata[ii:2,'Re']*Ygen[1:ii-1])/sum(Ygen[1:(ii-1)])
#     episimdata[ii, 'Lambda'] <- sum(episimdata[(ii-1):1,'I']*Ygen[1:(ii-1)])
#
#     pois_input <- sum(episimdata[(ii-1):1,'I']*episimdata[ii:2,'Re']*Ygen[1:(ii-1)])
#     #print(pois_input)
#     episimdata[ii,'I'] <- rpois(1, pois_input)
#
#     if (susceptibles == 1) {
#       episimdata[ii, 'S'] <- episimdata[(ii-1), 'S'] - episimdata[ii,'I']
#       if (episimdata[ii, 'S'] < 0) {
#         episimdata[ii, 'S'] = 0
#       }
#     }
#
#     if (delay == 1) {
#       pois_input_c <- sum(episimdata[ii:1,'I']* Ydel[1:ii])
#
#       episimdata[ii,'C'] <- rpois(1, pois_input_c)
#     }
#
#     if (ur == 1) {
#       if (delay == 1){
#         episimdata[ii,'C'] <- rbetabinom.ab(1, episimdata[ii,'C'], ur_beta_a, ur_beta_b)
#       } else {
#         episimdata[ii,'C'] <- rbetabinom.ab(1, episimdata[ii,'I'], ur_beta_a, ur_beta_b)
#       }
#     }
#
#     if (delay*ur != 1){
#       episimdata[ii,'C'] <- episimdata[ii,'I']
#     }
#
#
#   }
#   return (episimdata)
# }

for (jj in 1:sim_ens) {
  episim_data_ens[[jj]] <- Epi_run(episim_data_ens[[jj]], Epi_pars, Noise_pars, pathogen = 2, delay = 1, ur = 1, N = N)
}

#episim_data <- Epi_run(episim_data)

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
