#' A function for computing the variance estimation
#'
#' This function takes two vectors and returns an average of squared differences.
#'
#' @param a a vector of data
#' @param b a vector of data
#' @return The value of mean((a[i]-b[j])^2) taken over (i,j)'s.
#'


mean_ss <- function(a,b){
  # calculate the cross term difference squared
  t = 0
  for (i in 1:length(a)){
    for (j in 1:length(b)){
      t = t + (a[i] - b[j])^2
    }
  }
  return (t/length(a)/length(b))
}
