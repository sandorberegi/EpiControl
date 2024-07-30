# Clean the workspace and console
#closeAllConnections(); rm(list=ls())
#cat("\014"); graphics.off()

# Set working directory to source
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)
# Folder path for results
folres = paste0("./results/covid/")

# Main functions to run EpiFilter
#files.sources = list.files(path = "./epifilter_main")
#for (i in 1:length(files.sources)) {
#  source(paste0(c("./epifilter_main/", files.sources[i]), collapse = ''))
#}

library(haven)
data <- read_dta('covid_data/Data_from_OxCGRT.dta')

# Filter rows where country is "United Kingdom"
uk_data <- subset(data, country == "United Kingdom")

pop_data <- read_dta('covid_data/World_population_world_bank.dta')

uk_pop <- pop_data[251, 2]

# Display the filtered data
print(uk_data)

# Convert date column to Date type if necessary
uk_data$date_mdy <- as.Date(uk_data$date_mdy, format="%m/%d/%Y")

uk_data$daily_cases_100k <- round((uk_data$daily_cases_100k)*67215293/1e5, digits =  0)

uk_data$daily_cases_100k[is.na(uk_data$daily_cases_100k)] <- 0

first_nonzero_day <- which(uk_data$daily_cases_100k > 0)[1]

# Print the first nonzero day
print(first_nonzero_day)


real_days <- 370L

# Incidence and dates
Iday = uk_data$daily_cases_100k[first_nonzero_day:real_days]
dates  = uk_data$date_mdy[first_nonzero_day:real_days]
# Time series lengths
nday = real_days-first_nonzero_day+1; tday = 1:nday

gen_time <- 6.5
gen_time_var <- 2.1


# Approxumate serial interval distribution from Ferguson et al

#wdist = dgamma(tday, shape = 2.3669, scale = 2.7463)
wdist = dgamma(tday, gen_time/gen_time_var, 1/gen_time_var)

# Total infectiousness
Lday = rep(0, nday)
for(i in 2:nday){
  # Total infectiousness
  Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])
}

######################################################################
## EpiFilter: provides formally smoothed and exact estimates
# Method based on Bayesian recursive filtering and smoothing
######################################################################

# Setup grid and noise parameters
Rmin = 0.01; Rmax = 10; eta = 0.1

# Uniform prior over grid of size m
m = 200; pR0 = (1/m)*rep(1, m)
# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)

epiFilter <- function(Rgrid, m, eta, pR0, nday, Lday, Iday, a){

  # Probability vector for R and prior
  pR = matrix(0, nday, m); pRup = pR
  pR[1, ] = pR0; pRup[1, ] = pR0

  # Mean and median estimates
  Rmean = rep(0, nday); Rmed = Rmean
  # 50% and 95% (depends on a) confidence on R
  Rhat = matrix(0, 4, nday)

  # Initialise mean
  Rmean[1] = pR[1, ]%*%Rgrid
  # CDF of prior
  Rcdf0 = cumsum(pR0)
  # Initialise quartiles
  idm = which(Rcdf0 >= 0.5, 1); Rmed[1] = Rgrid[idm[1]]
  id1 = which(Rcdf0 >= a, 1); id2 = which(Rcdf0 >= 1-a, 1)
  id3 = which(Rcdf0 >= 0.25, 1); id4 = which(Rcdf0 >= 0.75, 1)
  Rhat[1, 1] = Rgrid[id1[1]]; Rhat[2, 1] = Rgrid[id2[1]]
  Rhat[3, 1] = Rgrid[id3[1]]; Rhat[4, 1] = Rgrid[id4[1]]

  # Precompute state distributions for R transitions
  pstate = matrix(0, m, m);
  for(j in 1:m){
    pstate[j, ] = dnorm(Rgrid[j], Rgrid, sqrt(Rgrid)*eta)
  }

  # Update prior to posterior sequentially
  for(i in 2:nday){
    # Compute mean from Poisson renewal (observation model)
    rate = Lday[i]*Rgrid
    # Probabilities of observations
    pI = dpois(Iday[i], rate)

    # State predictions for R
    pRup[i, ]  = pR[i-1, ]%*%pstate
    # Update to posterior over R
    pR[i, ] = pRup[i, ]*pI
    pR[i, ] = pR[i, ]/sum(pR[i, ])

    # Posterior mean and CDF
    Rmean[i] = pR[i, ]%*%Rgrid
    Rcdf = cumsum(pR[i, ])

    # Quantiles for estimates
    idm = which(Rcdf >= 0.5, 1); Rmed[i] = Rgrid[idm[1]]
    id1 = which(Rcdf >= a, 1); id2 = which(Rcdf >= 1-a, 1)
    id3 = which(Rcdf >= 0.25, 1); id4 = which(Rcdf >= 0.75, 1)
    Rhat[1, i] = Rgrid[id1[1]]; Rhat[2, i] = Rgrid[id2[1]]
    Rhat[3, i] = Rgrid[id3[1]]; Rhat[4, i] = Rgrid[id4[1]]
  }

  # Main outputs: estimates of R and states
  epiFilter = list(Rmed, Rhat, Rmean, pR, pRup, pstate)
}

# Filtered (causal) estimates as list [Rmed, Rhatci, Rmean, pR, pRup, pstate]
Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], 0.025)

recursPredict <- function(Rgrid, pR, Lday, Rmean, a){

  # Grid size and length of time series
  nday = nrow(pR); m = ncol(pR)
  # Test lengths of inputs
  if (length(Rgrid) != m | length(Lday) != nday){
    stop("Input vectors of incorrect dimension")
  }

  # Mean prediction: Lday[i] => Iday[i+1]
  pred = Lday*Rmean; pred = pred[1:length(pred)-1]

  # Discrete space of possible predictions
  Igrid = 0:80000; lenI = length(Igrid);

  # Check if close to upper bound
  if (any(pred > 0.9*max(Igrid))){
    stop("Epidemic size too large")
  }

  # Prediction cdf and quantiles (50% and 95%)
  Fpred = matrix(0, nday-1, lenI)
  predInt = matrix(0, 4, nday-1)

  pb <- txtProgressBar(min = 0, max = nday-1, initial = 0, style = 3, width = 30, char = "=")

  # At every time construct CDF of predictions
  for(i in 1:(nday-1)){
    # Compute rate from Poisson renewal
    rate = Lday[i]*Rgrid
    # Prob of any I marginalised over Rgrid
    pI = rep(0, lenI)

    # Probabilities of observations 1 day ahead
    for(j in 1:lenI){
      # Raw probabilities of Igrid
      pIset = dpois(Igrid[j], rate)
      # Normalised by probs of R
      pI[j] = sum(pIset*pR[i, ])
    }

    # Quantile predictions and CDF at i+1
    Fpred[i, ] = cumsum(pI)/sum(pI)
    id1 = which(Fpred[i, ] >= a); id2 = which(Fpred[i, ] >= 1-a)
    id3 = which(Fpred[i, ] >= 0.25); id4 = which(Fpred[i, ] >= 0.75)

    # Assign prediction results
    predInt[1, i] = Igrid[id1[1]]; predInt[2, i] = Igrid[id2[1]]
    predInt[3, i] = Igrid[id3[1]]; predInt[4, i] = Igrid[id4[1]]

    setTxtProgressBar(pb,i)
  }
  # Main outputs: mean and 95% predictions
  close(pb)
  recursPredict = list(pred, predInt)
}

# Causal predictions from filtered estimates [pred predci]
Ifilt = recursPredict(Rgrid, Rfilt[[4]], Lday[tday], Rfilt[[3]], 0.025)

library(ggplot2)

dfI <- data.frame(
  index = 1:length(Ifilt[[1]]),
  value1 = Ifilt[[1]]
)

ggplot(dfI) +
  geom_line(aes(x = index, y = value1, color = "Value 1")) +   # Line plot for value1
  geom_point(aes(x = index, y = value1, color = "Value 1"))# Add points for value

epiSmoother <- function(Rgrid, m, pR, pRup, nday, pstate, a){

  # Last smoothed distribution same as filtered
  qR = matrix(0, nday, m); qR[nday, ] = pR[nday, ]

  # Main smoothing equation iteratively computed
  for(i in seq(nday-1, 1)){
    # Remove zeros
    pRup[i+1, pRup[i+1, ] == 0] = 10^-8

    # Integral term in smoother
    integ = qR[i+1, ]/pRup[i+1, ]
    integ = integ%*%pstate

    # Smoothed posterior over Rgrid
    qR[i, ] = pR[i, ]*integ
    # Force a normalisation
    qR[i, ] = qR[i, ]/sum(qR[i, ]);
  }

  # Mean, median estimats of R
  Rmean = rep(0, nday); Rmed = Rmean
  # 50% and 95% (depends on a) confidence on R
  Rhat = matrix(0, 4, nday)

  # Compute at every time point
  for (i in 1:nday) {
    # Posterior mean and CDF
    Rmean[i] = qR[i, ]%*%Rgrid
    Rcdf = cumsum(qR[i, ])

    # Quantiles for estimates
    idm = which(Rcdf >= 0.5); Rmed[i] = Rgrid[idm[1]]
    id1 = which(Rcdf >= a, 1); id2 = which(Rcdf >= 1-a, 1)
    id3 = which(Rcdf >= 0.25, 1); id4 = which(Rcdf >= 0.75, 1)
    Rhat[1, i] = Rgrid[id1[1]]; Rhat[2, i] = Rgrid[id2[1]]
    Rhat[3, i] = Rgrid[id3[1]]; Rhat[4, i] = Rgrid[id4[1]]
  }

  # Main outputs: estimates of R and states
  epiSmoother = list(Rmed, Rhat, Rmean, qR)
}


# Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
Rsmooth = epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], nday, Rfilt[[6]], 0.025)
# Smoothed predictions from filtered estimates [pred predci]
Ismooth = recursPredict(Rgrid, Rsmooth[[4]], Lday[tday], Rsmooth[[3]], 0.025)

plotEpiFilter <- function(Rhat, Rhatci, Inexhat, Inexhatci, plotname, Iplt, folres, eta){

  # Check lengths
  if (length(Rhat) != length(Inexhat)){
    print(c('Rhat length', length(Rhat)))
    print(c('Ihat length', length(Inexhat)))
    stop('Inconsistent incidence and reprod. num vectors')
  }else{
    # Length of time (relative)
    tset = 1:length(Rhat)

    # Two panel plot of estimates and predictions
    pdf(file=paste0(folres, plotname, '.pdf'))
    par(mfrow=c(2,1))
    # Reprod. num estimates and confidence interval
    plot(tset, Rhat, type = 'l', bty = 'l', lwd = 2, col='blue',
         xlab = paste0("time (eta = ", eta, ")"), ylab = "reprod. number")
    polygon(c(tset, rev(tset)), c(Rhatci[1,], rev(Rhatci[2,])),
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Rhatci[3,], rev(Rhatci[4,])),
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    lines(tset, rep(1, length(tset)), lwd = 2, col = 'black', lty = 'dashed')

    # Incidence predictions and confidence interval
    plot(tset, Inexhat, type = 'l', bty = 'l', lwd = 2, col='blue',
         xlab = paste0("time (eta = ", eta, ")"), ylab = "incidence", ylim = c(0, max(Iplt)+30))
    polygon(c(tset, rev(tset)), c(Inexhatci[1,], rev(Inexhatci[2,])),
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Inexhatci[3,], rev(Inexhatci[4,])),
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    points(tset, Iplt, pch = 19, col = 'gray')
    dev.off()

  }
}

# Plot estimates and predictions from filtering
plotEpiFilter(Rfilt[[3]][2:nday], Rfilt[[2]][, 2:nday], Ifilt[[1]], Ifilt[[2]],
              'EpiFilter_UK', Iday[2:nday], folres, eta)

# Plot estimates and predictions from smoothing
plotEpiFilter(Rsmooth[[3]][2:nday], Rsmooth[[2]][, 2:nday], Ismooth[[1]], Ismooth[[2]],
              'EpiSmooth_UK', Iday[2:nday], folres, eta)


Rfilt_estimates = Rfilt[[3]]
Rsmooth_estimates = Rsmooth[[3]]
library(ggplot2)
# Create a data frame for plotting
data = data.frame(
  Time = rep(tday, 2),
  ReproductionNumber = c(Rfilt_estimates, Rsmooth_estimates),
  Type = factor(rep(c('Filtered', 'Smoothed'), each = length(tday)))
)

# Create the plot using ggplot2
ggplot(data, aes(x = Time, y = ReproductionNumber, color = Type)) +
  geom_line() +
  labs(title = 'Reproduction Number Estimates', x = 'Time (days)', y = 'Reproduction Number') +
  theme_minimal() +
  scale_color_manual(values = c('Filtered' = 'blue', 'Smoothed' = 'red'))

last_Rs <- tail(Rsmooth[[3]], 20)

Rest_const <- mean(last_Rs)

library(fitdistrplus)

fit <- fitdist(last_Rs, "gamma")
fit$estimate

dfR <- data.frame(value = last_Rs)

ggplot(dfR, aes(x = value)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.2, fill = "blue", color = "black", alpha = 0.7) +
  stat_function(fun = dgamma, args = list(shape = fit$estimate['shape'], rate = fit$estimate['rate']),
                color = "red", size = 1) +
  labs(title = "Histogram and Fitted Gamma Distribution",
       x = "Value",
       y = "Density") +
  theme_minimal()


### Run my code from here

library(VGAM)
#library(future.apply)
#library(foreach)
library(parallel)
library(pbapply)

cores=detectCores()-1
cl <- makeCluster(cores)

Epi_pars <- data.frame (
  Pathogen = c("COVID-19", "Ebola"),
  R0 = c(Rest_const, 2.5),
  gen_time = c(6.5, 15.0),
  gen_time_var = c(2.1, 2.1),
  R0_shape = c(fit$estimate["shape"], 1.0),
  R0_rate = c(fit$estimate["rate"], 1.0)
)

#real_days <- 258L #336
real_days <- 370L

I0 <- 10L #initial no. of infections
ndays <- real_days+41L*7L#epidemic length
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
column_names <- c("days", "sim_id", "I", "Lambda", "C", "Lambda_C", "S", "Re", "Rew", "Rest", "R0est", "policy", "R_coeff", "Real_C")

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

episim_data['sim_id'] <- rep(1, ndays)
episim_data[1:real_days,] <- c(1, 1, I0, I0, Noise_pars['ur_mean']*I0, Noise_pars['ur_mean']*I0, N-I0, Epi_pars[1,'R0'], Epi_pars[1,'R0'], 1, 1, 1, 1, 0)
episim_data['days'] <- 1:ndays
episim_data['policy'] <- rep(1, ndays)

episim_data[1:real_days,"C"] <- uk_data$daily_cases_100k[1:real_days]
episim_data[1:real_days,"I"] <- uk_data$daily_cases_100k[1:real_days]
episim_data[1:nrow(episim_data),"Real_C"] <- uk_data$daily_cases_100k[1:nrow(episim_data)]

# get infectiousness and estimated R-s

gen_time <- 6.5
gen_time_var <- 2.1

R_est_wind <- 5  # Define your window size

Ygen <- dgamma(1:nrow(uk_data), gen_time/gen_time_var, 1/gen_time_var)
Ygen <- Ygen/sum(Ygen)

# Estimate R
for (ii in 1:real_days+1) {
  if (ii-1 < R_est_wind) {
    episim_data[ii, 'Rest'] <- mean(episim_data[1:(ii-1), 'C']) / mean(episim_data[1:(ii-1), 'Lambda_C'])
    R_coeff_tmp <- sum(Ygen[1:(ii-1)] * episim_data[(ii-1):1, 'R_coeff']) / sum(Ygen[1:(ii-1)])
  } else {
    if ( mean(episim_data[(ii-R_est_wind):(ii-1), 'Lambda_C']) == 0){
      episim_data[ii, 'Rest'] <- 0
      R_coeff_tmp <- 1
    } else {
      episim_data[ii, 'Rest'] <- mean(episim_data[(ii-R_est_wind):(ii-1), 'C']) / mean(episim_data[(ii-R_est_wind):(ii-1), 'Lambda_C'])
      R_coeff_tmp <- sum(Ygen[1:(ii-1)] * episim_data[(ii-1):1, 'R_coeff']) / sum(Ygen[1:(ii-1)])
    }
  }
  episim_data[ii, 'R0est'] <- episim_data[ii, 'Rest'] / R_coeff_tmp
  episim_data[ii, 'Lambda_C'] <- sum(episim_data[(ii-1):1,'C']*Ygen[1:(ii-1)])
  episim_data[ii, 'Lambda'] <- sum(episim_data[(ii-1):1,'I']*Ygen[1:(ii-1)])
}

#Epi_pars[1,"R0"] <- episim_data[real_days, 'R0est']

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
clusterEvalQ(cl, Epi_pars)
clusterEvalQ(cl, {
  library(VGAM)
  library(EpiControl)
})

#for (ii in 1:sim_ens) {
#  episim_data_ens[[jj]] <- Epi_MPC_run(episim_data_ens[[jj]], Epi_pars, Noise_pars, Action_space, pred_days = pred_days, n_ens = n_ens, start_day = real_days, ndays = nrow(episim_data), R_est_wind = R_est_wind, pathogen = 1, susceptibles = 0, delay = 0, ur = 0, N = N)
#}




results <- pblapply(1:sim_ens, function(jj) {
  episim_data_ens[[jj]] <- Epi_MPC_run_Rdistr(episim_data_ens[[jj]], Epi_pars, Noise_pars, Action_space, pred_days = pred_days, n_ens = n_ens, start_day = real_days, ndays = nrow(episim_data), R_est_wind = R_est_wind, pathogen = 1, susceptibles = 0, delay = 0, ur = 0, R_uncert = 1, c_uncert = 1, N = N)
}, cl=cl)

episim_data_ens <- results
stopCluster(cl)



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

filtered_data1 <- subset(combined_data, days > real_days)
filtered_data2 <- subset(combined_data, days <= real_days)

ggplot() +
  geom_line(data = subset(filtered_data1, sim_id != 20), aes(x = days, y = Rest, color = as.factor(sim_id)), alpha = 0.3) +
  geom_line(data = subset(filtered_data1, sim_id == 20), aes(x = days, y = Rest), color = "red", alpha = 1.0) +
  geom_line(data = subset(filtered_data2, sim_id != 20), aes(x = days, y = Rest, color = as.factor(sim_id)), alpha = 0.3) +
  geom_line(data = subset(filtered_data2, sim_id == 20), aes(x = days, y = Rest), color = "red", alpha = 1.0) +
  labs(x = "Days", y = "R0 estimate") +
  scale_color_manual(values = cls) +
  guides(color = FALSE)

ggplot() +
  geom_line(data = subset(filtered_data1, sim_id != 20), aes(x = days, y = Re, color = as.factor(sim_id)), alpha = 0.3) +
  geom_line(data = subset(filtered_data1, sim_id == 20), aes(x = days, y = Re), color = "red", alpha = 1.0) +
  geom_line(data = subset(filtered_data2, sim_id != 20), aes(x = days, y = Rest, color = as.factor(sim_id)), alpha = 0.3) +
  geom_line(data = subset(filtered_data2, sim_id == 20), aes(x = days, y = Rest), color = "red", alpha = 1.0) +
  geom_vline(xintercept = real_days, linetype = "dashed", color = "black", size=0.25) +
  labs(x = "Days", y = "R") + xlim(50,650) + ylim(0,8) +
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
ggplot(combined_data %>% filter(sim_id == 20)) +
  geom_line(data = subset(combined_data, sim_id != 20), aes(x = days, y = C, color = as.factor(sim_id)), alpha = 0.1) +
  geom_line(aes(x = days, y = C, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0) +
  geom_line(aes(x = days, y = I, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0, size=0.25) +
  geom_line(aes(x = days, y = Real_C, color = "blue", group = 1), alpha = 1.0, size=0.25) +
  geom_hline(yintercept = C_target, linetype = "dashed", color = "blue", size=0.25) +
  geom_vline(xintercept = real_days, linetype = "dashed", color = "black", size=0.25) +
  labs(x = "Days", y = "Reported cases and true infections", color = "Policy") +
  scale_color_manual(values = c("No intervention" = "chartreuse3", "Social distancing" = "darkorchid1", "Lockdown" = "red")) +
  guides(color = guide_legend(title = "Policy"))

