#' A function for covariance estimation (one component of the 3rd term)
#'
#' This function takes a few arguments and compute cone compunent of the
#' third term of the covariance bound estimator
#'
#' @param a a vector of data, gamma' of a treatment group l
#' @param b a vector of data, gamma' of another treatment group m
#' @param pil a real number, the propensity score of treatment l
#' @param pim a real number, the propensity score of treatment m
#' @param nb an integer, the block size
#' @return the value of mean(a[i]b[j]) / pil / pim
#'


thd_term_1 <- function(a,b,pil,pim,nb){
  # calculate the third term in the covariance est 3
  t = 0
  nbl = length(a)
  nbm = length(b)
  for (i in 1:nbl){
    for (j in 1:nbm)
      t = t + a[i] * b[j]
  }
  t = nb * (nb-1) / nbl / nbm * (nbl / nb) * nbm / (nb-1) / pil / pim * t
  return (t)
}
