#' Hajék estimator with sandwich-type SEs
#'
#' This function takes a dataframe and returns the Hajék estimator with
#' 3 or 1 variance estimator and 4 or 3 covariance estimator depending on the
#' block and treatment group sizes
#'
#' @param data a dataframe of weights, outcomes, treatment assignments, and block indicators
#' @param weight a character, the name of the weight column
#' @param y a character, the name of the outcome column
#' @param z a character, the name of the treatment indicator column
#' @param block a character, the name of the block indicator column
#' @return a list of results. $theta: a vector of Hajék ratios; $varest: a matrix of varaince estimates, $covest: an array of covariance bound estimates.
#'
hajek_est <- function(data, weight, y, z, block){
  ws <- data[[weight]]  # weights

  yobs <- data[[y]]  # observed ys
  zobs <- data[[z]]  # observed zs (random assignment)
  nbs <- table(data[[block]])  # block sizes
  # names(nbs) <- NULL

  nbk <- table(data.frame(data[[block]], data[[z]]))
  # nbk, number of units in each block and treatment, B by K
  nbk_all <- nbk[data[[block]], ]
  # nbk, number of units in each block and treatment
  # corresponding to each unit, n by K

  B <- length(nbs)  # number of blocks
  K <- ncol(nbk)  # number of treatments
  mu <- sum(ws) / B  # mu_U

  pi <- t(t(nbk) / t(nbs)[1,])
  # pi_b,k propensity score of each block and treatment, B by K
  pi_all <- pi[data[[block]], ]
  # pi_b,k propensity score of each block and treatment
  # corresponding to each unit, n by K

  gammas <- cbind(nbk_all[,1] * ws, nbk_all[,2] * ws, nbk_all[,3] * ws)
  # gamma (or xps and yps) for variance estimation, n by K
  gamsbk <- matrix(nrow=B, ncol=K)  # s^2_b,j, sample variance
  sigmab <- matrix(nrow=B, ncol=K)  # hat sigma b,j

  nu1 <- matrix(nrow=B, ncol=K-1)  # nu1_b,0k
  nu2 <- matrix(nrow=B, ncol=K-1)  # nu2_b,0k
  nu3 <- matrix(nrow=B, ncol=K-1)  # nu3_b,0k

  gammas_cov <- cbind(ws, ws, ws)  # gamma' for covariance estimation, n by K
  bd2 <- array(dim=c((K-1)*(K-2)/2, B, 2))  # bound 2 0,k,j
  bd3 <- array(dim=c((K-1)*(K-2)/2, B, 1))  # bound 3 0,k,j
  bd4 <- array(dim=c((K-1)*(K-2)/2, B, 2))  # bound 4 0,k,j

  theta <- vector(length=K)  # Hajek estimators
  varest <- matrix(nrow=3, ncol=K-1)  # variance estimators
  covest <- array(dim=c(K-1, K-1, 4))  # covariance estimators

  for (k in 1:K){
    indk <- (zobs == k)
    tk <- sum(ws[indk] * yobs[indk] / pi_all[indk, k]) /
      sum(ws[indk] / pi_all[indk, k])  # ratio estimates
    theta[k] <- tk

    gammas[indk,k] <- gammas[indk,k] / pi_all[indk, k] * (yobs[indk] - theta[k])
    gammas_cov[indk,k] <- gammas_cov[indk,k] * (yobs[indk] - theta[k])

    for (b in 1:B){
      in_b <- (sum(nbs[1:b])-nbs[b]+1):(sum(nbs[1:b]))  # indices of units in block b
      indbk <- (zobs[in_b] == k)  # indices of units in block b, treatment k

      if (sum(indbk) == 1)
        gamsbk[b,k] = 0
      else
        gamsbk[b,k] = sd(gammas[in_b,k][indbk])^2

      sigmab[b,k] <- gamsbk[b,k] * (nbs[b] - nbs[b]*pi[b,k]) / nbs[b]^2 / pi[b,k]

      # variance estimators by block
      if (k > 1){
        indb0 <- (zobs[in_b] == 1)

        nu3[b,k-1] <- mean_ss(gammas[in_b,k][indbk],
                              gammas[in_b,1][indb0]) -
          (nbk[b,1]-1) / nbk[b,1] * gamsbk[b,1] -
          (nbk[b,k]-1) / nbk[b,k] * gamsbk[b,k]

        if (sigmab[b,1] == 0 | sigmab[b,k] == 0){
          # in this case, nu1 and nu2 are invalid, fill in with nu3 instead
          nu1[b,k-1] <- nu3[b,k-1]
          nu2[b,k-1] <- nu3[b,k-1]
        }
        else{
          nu1[b,k-1] <- nbs[b] / (nbs[b] - nbk[b,1]) * sigmab[b,1] +
            nbs[b] / (nbs[b] - nbk[b,k]) * sigmab[b,k]

          nu2[b,k-1] <- (1 + nbk[b,k] / (nbs[b] - nbk[b,1])) * sigmab[b,1] +
            (1 + nbk[b,1] / (nbs[b] - nbk[b,k])) * sigmab[b,k]
        }
      }

      # covariance estimators by block
      if (k > 2){
        for (j in 2:(k-1)){
          indbj <- (zobs[in_b] == j)
          m <- sum(0:(k-3)) + j - 1  # index for the bound of cov_k,j

          bd2[m,b,1] <- -
            mean_ss(gammas_cov[in_b,k][indbk]/pi[b,k],
                    gammas_cov[in_b,1][indb0]/pi[b,1]) /
            (2 + 2 / (nbs[b] - 1)) -
            mean_ss(gammas_cov[in_b,j][indbj]/pi[b,j],
                    gammas_cov[in_b,1][indb0]/pi[b,1]) /
            (2 + 2 / (nbs[b] - 1)) -
            mean_ss(gammas_cov[in_b,k][indbk]/pi[b,k],
                    gammas_cov[in_b,j][indbj]/pi[b,j]) /
            (2 - 2 / (nbs[b] - 1))

          bd2[m,b,2] <-
            mean_ss(gammas_cov[in_b,k][indbk]/pi[b,k],
                    gammas_cov[in_b,1][indb0]/pi[b,1]) /
            (2 - 2 / (nbs[b] - 1)) +
            mean_ss(gammas_cov[in_b,j][indbj]/pi[b,j],
                    gammas_cov[in_b,1][indb0]/pi[b,1]) /
            (2 - 2 / (nbs[b] - 1)) +
            mean_ss(gammas_cov[in_b,k][indbk]/pi[b,k],
                    gammas_cov[in_b,j][indbj]/pi[b,j]) /
            (2 + 2 / (nbs[b] - 1))

          if (sum(sigmab[b,] == 0) > 0)
            bd3[m,b,1] <- 0
          else
            bd3[m,b,1] <- sqrt(sigmab[b,1]*sigmab[b,j]) +
            sqrt(sigmab[b,1]*sigmab[b,k]) + sqrt(sigmab[b,j]*sigmab[b,k])

          bd4[m,b,1] =
            scd_term(gammas_cov[in_b,1][indb0], nbs[b]) +
            scd_term(gammas_cov[in_b,j][indbj], nbs[b]) +
            scd_term(gammas_cov[in_b,k][indbk], nbs[b])
          bd4[m,b,2] <- thd_term(
            gammas_cov[in_b,1][indb0],
            gammas_cov[in_b,j][indbj],
            gammas_cov[in_b,k][indbk],
            pi[b,1],
            pi[b,j],
            pi[b,k],
            nbs[b]
          )
        }
      }
    }

    # variance estimators
    if (k > 1){
      varest[1,k-1] <- mean(nu1[,k-1]) / mu^2
      varest[2,k-1] <- mean(nu2[,k-1]) / mu^2
      varest[3,k-1] <- mean(nu3[,k-1]) / mu^2
    }
    # covariance estimators
    if (k > 2){
      for (j in 2:(k-1)){
        m <- sum(0:(k-3)) + j - 1  # index for the bound of cov_k,j

        if (sum(sigmab[,1] == 0) > 0)  # estimate sigma_0^2
          sigmab0 <- mean(gammas[,1]^2)
        else
          sigmab0 <- mean(sigmab[,1])

        covest[j-1, k-1, 1] <- - sqrt(varest[1,j-1] * varest[1,k-1])
        covest[k-1, j-1, 1] <- - covest[j-1, k-1, 1]

        covest[j-1, k-1, 2] <- (sigmab0 + mean(bd2[m,,1])) / mu^2
        covest[k-1, j-1, 2] <- (sigmab0 + mean(bd2[m,,2])) / mu^2

        covest[j-1, k-1, 3] <- (sigmab0 - mean(bd3[m,,1])) / mu^2
        covest[k-1, j-1, 3] <- (sigmab0 + mean(bd3[m,,1])) / mu^2

        covest[j-1, k-1, 4] <- (sigmab0 - mean(bd4[m,,1]) + mean(bd4[m,,2])) / mu^2
        covest[k-1, j-1, 4] <- (sigmab0 + mean(bd4[m,,1]) + mean(bd4[m,,2])) / mu^2
      }
    }

  }
  if (sum(bd3 == 0) > 0)
    return(list(theta=theta, varest=varest/B, covest=covest[,,c(1,2,4)]/B))
  return(list(theta=theta, varest=varest/B, covest=covest/B))
}
