softmax <- function(x) {
  exp_x <- exp(x)  # Exponentiate the input vector
  return(exp_x / sum(exp_x))  # Normalize by dividing by the sum of all exponentiated values
}
