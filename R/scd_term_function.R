#' A function for covariance estimation (2nd term)
#'
#' This function takes a vector and an integer and returns the second
#' term of the covariance bound estimator
#'
#' @param a a vector of data, gamma's of a block
#' @param nb an integer, the block size
#' @return the value of nb * (mean(a[i]^2) + mean(a[i]a[j])/(nb-1))
#'

scd_term <- function(a,nb){
  # calculate the second term in the covariance est 3
  t = 0
  nbl = length(a)
  for (i in 1:nbl){
    t = t + a[i]^2
    for (j in 1:nbl){
      if (j != i){
        t = t + (nb-1) / (nbl-1) * a[i]*a[j]
      }
    }
  }
  return (t * nb / nbl)
}
