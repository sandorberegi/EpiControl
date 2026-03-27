#' Runs MPC while stepping a custom compartmental model
#' Requires:
#'   - model: built via build_model()
#'   - step_model(), state(), flow_SIRDS() that attaches events via attr(delta, "events")
#'   - episimdata: has columns for compartments and events
#'   - R_estimator(): function to estimate R
#'   - action space
#'   - updatepars: the function defining how interverntions affect transmission
#'   - disease parameters
#'
#' @export

Epi_MPC_run_comp_model <- function(episimdata,
                              episettings, epi_par, noise_par, actions,
                              pred_days, n_ens = 100,
                              start_day = 1,
                              ndays = nrow(episimdata),
                              R_est_wind = 5, pathogen = 1, susceptibles = 0,
                              delay = 0, ur = 0, r_dir = 1, N = sim_settings$N) {

  model <- episettings$model

  R_estimator <- episettings$R_estimator
  sim_settings <- episettings$sim_settings

  rf           <- sim_settings$rf
  r_trans_steep<- sim_settings$r_trans_steep
  t0           <- sim_settings$t0

  # Epidemiological parameters for the renewal model (for projections)
  R0        <- epi_par[pathogen, "R0"]
  gen_time  <- epi_par[pathogen, "gen_time"]
  gen_time_var <- epi_par[pathogen, "gen_time_var"]

  # Generation interval kernel (Gamma distribution)
  Ygen <- dgamma(1:ndays, shape = gen_time / gen_time_var, rate = 1 / gen_time_var)
  Ygen <- Ygen / sum(Ygen)

  number_of_actions <- nrow(actions)

  # Reporting delay (optional)
  if (delay == 1) {
    repd_mean <- noise_par[1, 'repd_mean']
    del_disp  <- noise_par[1, 'del_disp']
    Ydel <- dgamma(1:ndays, shape = del_disp, rate = del_disp / repd_mean)
    Ydel <- Ydel / sum(Ydel)
  }

  # Under-reporting (optional)
  if (ur == 1) {
    ur_mean   <- noise_par[1, 'ur_mean']
    ur_beta_a <- noise_par[1, 'ur_beta_a']
    ur_beta_b <- (1 - ur_mean) / ur_mean * ur_beta_a
  }

  for (ii in (start_day + 1):ndays) {

    ## --- 1) Estimate reproduction number from data up to previous day
    R_est_res <- R_estimator(episimdata, Ygen, ii, R_est_wind = R_est_wind, r_dir = r_dir)
    episimdata[ii, 'Rest']    <- R_est_res$R_est
    R_coeff_tmp               <- R_est_res$R_coeff_tmp
    episimdata[ii, 'R0est']   <- episimdata[ii, 'Rest'] / R_coeff_tmp

    ## --- 2) Choose policy every rf days by expected reward
    if (ii %% rf == 0L) {
      Rewards <- numeric(number_of_actions)
      for (jj in 1:number_of_actions) {
        Reward_ens <- replicate(n_ens, Epi_pred_wd(
          episimdata, episettings, epi_par, noise_par, actions,
          pathogen, pred_days, r_dir, ii, jj, N, ndays = ndays
        ))
        Rewards[jj] <- mean(Reward_ens)
      }
      episimdata[ii, 'policy'] <- which.max(Rewards)
    } else {
      episimdata[ii, 'policy'] <- episimdata[ii - 1, 'policy']
    }

    ## --- 3) Compute effective R, and total infectiousness
    Rcoeff <- actions[episimdata[ii, 'policy'], 'R_coeff']
    episimdata[ii, 'R_coeff'] <- Rcoeff

    Re_act <- Rcoeff * R0
    episimdata[ii, 'Re']  <- Re_act

    episimdata[ii, 'Rew'] <- sum(episimdata[ii:2, 'Re'] * Ygen[1:(ii - 1)]) / sum(Ygen[1:(ii - 1)])

    # incidence-based renewal predictors
    episimdata[ii, 'Lambda']    <- sum(episimdata[(ii - 1):1, 'I'] * Ygen[1:(ii - 1)])
    episimdata[ii, 'Lambda_C']  <- sum(episimdata[(ii - 1):1, 'C'] * Ygen[1:(ii - 1)], na.rm = TRUE)

    if ((r_dir == 2) && (ii > rf)) {
      Rdir <- logistic_function((ii %% rf),
                                episimdata[(ii - (ii %% rf) - 1), 'Re'],
                                episimdata[ii, 'Re'],
                                r_trans_steep, t0)
      episimdata[ii, 'Rew'] <- Rdir
      Re_use <- Rdir
    } else {
      Re_use <- Re_act
    }

    pars_today <- episettings$model_par

    # Create a generic context object with anything computed in the step that may be used in updatepars().
    ctx <- list(
      ii = ii,
      t = ii - 1,
      dt = 1,
      policy = episimdata[ii, "policy"],
      Rcoeff = Rcoeff,
      Re_act = Re_act,
      Re_use = Re_use,
      Rest = episimdata[ii, "Rest"],
      R0est = episimdata[ii, "R0est"],
      Lambda = episimdata[ii, "Lambda"],
      Lambda_C = episimdata[ii, "Lambda_C"],
      episim_row_prev = if (ii > 1) episimdata[ii - 1, ] else NULL,
      sim_settings = sim_settings,
      episettings = episettings,
      epi_par = epi_par,
      noise_par = noise_par,
      actions = actions,
      pathogen = pathogen
    )

    # Apply user-defined intervention
    if (!is.null(episettings$updatepars)) {
      pars_today <- episettings$updatepars(pars_today, ctx)
      if (!is.list(pars_today)) stop("episettings$updatepars must return a list.")
    }

    step_res <- step_model(model, t = ii - 1, dt = 1, pars = pars_today,
                           event_names = episettings$event_names)


    model <- step_res$model

    st <- step_res$states
    episimdata[ii, names(st)] <- as.numeric(st)

    episimdata[ii, event_names] <- as.numeric(step_res$events)

    episimdata[ii, "I"] <- episimdata[ii, episettings$incidence_flow]

    # Observation model for reported cases I -> C (delay / under-reporting)
    if (delay == 1) {
      pois_input_c <- sum(episimdata[ii:1, 'I'] * Ydel[1:ii])
      episimdata[ii, 'C'] <- rpois(1, pois_input_c)
    }

    if (ur == 1) {
      if (delay == 1) {
        episimdata[ii, 'C'] <- VGAM::rbetabinom.ab(1, episimdata[ii, 'C'], ur_beta_a, ur_beta_b)
      } else {
        episimdata[ii, 'C'] <- VGAM::rbetabinom.ab(1, episimdata[ii, 'I'], ur_beta_a, ur_beta_b)
      }
    }

    if (delay + ur == 0) {
      episimdata[ii, 'C'] <- episimdata[ii, 'I']
    }
  }

  return(episimdata)
}
