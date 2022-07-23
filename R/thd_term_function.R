#' A function for covariance estimation (3rd term)
#'
#' This function takes a few arguments and compute the third term
#' of the covariance bound estimator.
#'
#' @param a a vector of data, gamma' of a treatment group 1
#' @param b a vector of data, gamma' of another treatment group 2
#' @param c a vector of data, gamma' of another treatment group 3
#' @param pi1 a real number, the propensity score of treatment 1
#' @param pi2 a real number, the propensity score of treatment 2
#' @param pi3 a real number, the propensity score of treatment 3
#' @param nb an integer, the block size
#' @return the value of - mean(a[i]b[j])/pi1/pi2 - mean(a[i]c[j])/pi1/pi3 + mean(b[i]c[j])/pi2/pi3

thd_term <- function(a,b,c,pi1,pi2,pi3,nb){
  # calculate the third term in the covariance est 3
  t = - thd_term_1(a,b,pi1,pi2,nb) - thd_term_1(a,c,pi1,pi3,nb) +
    thd_term_1(b,c,pi2,pi3,nb)
  return (t)
}
