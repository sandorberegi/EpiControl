library(EpiControl)
library(VGAM)
library(parallel)
library(pbapply)
library(zoo)  # For rolling sum operations
library(dplyr)
library(ggplot2)

Epi_pars <- data.frame(
  Pathogen = c("COVID-19", "Ebola"),
  R0 = c(2.2, 2.5),
  gen_time = c(6.5, 15.0),
  gen_time_var = c(2.1, 2.1),
  CFR = c(0.0132, 0.5),
  mortality_mean = c(10.0, 10.0),
  mortality_var = c(1.1, 1.1)
)

Noise_pars <- data.frame(
  repd_mean = 10.5,  # Reporting delay mean
  del_disp = 5.0,    # Reporting delay variance
  ur_mean = 0.3,     # Under-reporting mean
  ur_beta_a = 50.0   # Beta distribution alpha for under-reporting
)

Action_space <- data.frame(
  NPI = c("No restrictions", "Lockdown"),
  R_coeff = c(1.0, 0.3),
  R_beta_a = c(0.0, 5.0),
  cost_of_NPI = c(0.0, 0.15)
)

C_target <- 5000
D_target <- 12
r_trans_len <- 7

sim_settings <- list(
  ndays = 10 * 7, #simulation length
  start_day = 1,
  N = 1e7, # population size
  I0 = 10, # initial infections
  C_target = C_target, #target cases
  C_target_pen = C_target*1.5, #overshoot penalty threshold
  R_target = 1.0,
  D_target = D_target, #one way to get peaks at 400 is to increase this to 15
  D_target_pen = 50, #max death
  alpha = 0*1.3/C_target, #~proportional gain (regulates error in cases) covid
  #alpha = 3.25/C_target #~proportional gain (regulates error in cases) ebola
  alpha_d = 1.3/D_target,
  ovp = 0*5.0, #overshoot penalty
  dovp = 10.0, #death overshoot penalty
  gamma = 0.95, #discounting factor
  n_ens = 20L, #MC assembly size for 4
  sim_ens = 10L, #assembly size for full simulation
  rf = 21L, #days 14
  R_est_wind = 5L, #rf-2 #window for R estimation
  susceptibles = 0,
  pred_days = 28L,
  r_trans_steep = 1.5,  # Growth rate
  r_trans_len = r_trans_len,  # Number of days for the transition
  t0 = r_trans_len  / 2, # Midpoint of the transition
  pathogen = 1,
  susceptibles = 0,
  delay = 1,
  ur = 1,
  r_dir = 2,
  LD_on = 14, #on threshold
  LD_off = 7, #off threshold
  v_max_rate = 0.8,
  vac_scale = 100,
  vac_start = 370,
  delta_scale = 40,
  delta_start = 550,
  delta_multiplier = 1.75,
  v_protection_delta = (58+85)/200,
  v_protection_alpha = 0.83
)

episettings <- list(
  sim_function = Epi_MPC_run_V,
  reward_function = reward_fun,
  R_estimator = R_estim,
  noise_par = Noise_pars,
  epi_par = Epi_pars,
  actions = Action_space,
  sim_settings = sim_settings,
  parallel = FALSE
)

column_names <- c("days", "sim_id", "I", "Lambda", "C", "Lambda_C", "S", "Deaths", "Re", "Rew", "Rest", "R0est", "policy", "R_coeff")
episim_data <- data.frame(matrix(0, nrow = (sim_settings$ndays), ncol = length(column_names)))
colnames(episim_data) <- column_names

episim_data$policy <- rep(1, sim_settings$ndays)
episim_data$days <- 1:(sim_settings$ndays)
episim_data[1, ] <- c(1, 1, sim_settings$I0, sim_settings$I0, Noise_pars$ur_mean * sim_settings$I0, Noise_pars$ur_mean * sim_settings$I0, sim_settings$N - sim_settings$I0, 0, Epi_pars[1, "R0"], Epi_pars[1, "R0"], 1, 1, 1, 1)

episim_data_ens <- replicate(sim_settings$sim_ens, episim_data, simplify = FALSE)

for (ii in 1:sim_settings$sim_ens) {
  episim_data_ens[[ii]]$sim_id <- rep(ii, sim_settings$ndays)
}


cores=detectCores()-1
cl <- makeCluster(cores)

clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
clusterEvalQ(cl, library(EpiControl))
clusterEvalQ(cl, library(VGAM))

episettings$cl <- cl

results <- epicontrol(episim_data_ens, episettings)

stopCluster(cl)


