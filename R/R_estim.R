R_estim <- function(episimdata, Ygen, ii, R_est_wind = 5, r_dir = 0){
  if (ii-1 < R_est_wind) {
    R_est <- mean(episimdata[1:(ii-1), 'C'])/mean(episimdata[1:(ii-1), 'Lambda_C'])
    R_coeff_tmp <- sum(Ygen[1:(ii-1)] * episimdata[(ii-1):1, 'R_coeff'])/sum(Ygen[1:(ii-1)])
  } else {

    R_est <- mean(episimdata[(ii-R_est_wind):(ii-1), 'C'])/mean(episimdata[(ii-R_est_wind):(ii-1), 'Lambda_C'])

    if (r_dir == 1){
      R_coeff_tmp <-  mean(episimdata[(ii-R_est_wind):(ii-1), 'R_coeff'])
    } else {
      R_coeff_tmp <- sum(Ygen[1:(ii-R_est_wind)] * episimdata[(ii-R_est_wind):1, 'R_coeff'])/sum(Ygen[1:(ii-R_est_wind)])
    }
  }

  return(list(R_est = R_est, R_coeff_tmp = R_coeff_tmp))
}
