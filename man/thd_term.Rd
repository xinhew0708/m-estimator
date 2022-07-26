% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thd_term_function.R
\name{thd_term}
\alias{thd_term}
\title{A function for covariance estimation (3rd term)}
\usage{
thd_term(a, b, c, pi1, pi2, pi3, nb)
}
\arguments{
\item{a}{a vector of data, gamma' of a treatment group 1}

\item{b}{a vector of data, gamma' of another treatment group 2}

\item{c}{a vector of data, gamma' of another treatment group 3}

\item{pi1}{a real number, the propensity score of treatment 1}

\item{pi2}{a real number, the propensity score of treatment 2}

\item{pi3}{a real number, the propensity score of treatment 3}

\item{nb}{an integer, the block size}
}
\value{
the value of - mean(a[i]b[j])/pi1/pi2 - mean(a[i]c[j])/pi1/pi3 + mean(b[i]c[j])/pi2/pi3
}
\description{
This function takes a few arguments and compute the third term
of the covariance bound estimator.
}
