##' Generate Survival Data
##' @rdname genData
##' @usage genData(lambda, beta, n, tmax, epsilon)
##' 
##' @param lambda baseline hazard rate for exponential baseline
##' @param beta vector of coefficients
##' @param n number of observations
##' @param tmax maximum follow-up time
##' @param epsilon controls the censoring rate
##' 
##' @return a dataframe of five columns:
##' * \code{survtime} - survival times
##' * \code{status} - final statuses
##' * \code{age} - center and scaled age
##' * \code{Sex} - gender indicator
##' * \code{Black} - race indicator
##' @importFrom stats quantile rexp
##' @importFrom utils tail 
##' @export

genData <- function(lambda, beta, n, tmax, epsilon) {
  beta <- do.call(matrix, args = list(beta, length(beta)))
  rtmax <- runif(n, 0, tmax)
  vepsilon <- rbinom(n, 1, epsilon)
  rtmax <- tmax * vepsilon + (1 - vepsilon) * rtmax
  rtmax <- pmin(rtmax, tmax)
  p <- length(beta)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  x3 <- rbinom(n, 1, 0.1)
  x <- cbind(x1, x2, x3)
  hazard <- lambda * exp(x %*% beta)
  tfail <- rexp(n, rate = 1)
  tfail <- tfail / hazard
  status <- (tfail <= rtmax) + 0
  time <- tfail * status + rtmax * (1 - status)
  return(data.frame(
    survtime = time,
    status = status,
    age = x1,
    Sex = x2,
    Black = x3
  ))
}

##' Generate Survival Data
##' @rdname genDataWeibull
##' @usage genDataWeibull(weibA, weibB, beta, n, tmax, epsilon)
##' 
##' @param weibA Weibull shape parameter
##' @param weibB Weibull rate parameter
##' @param beta vector of coefficients
##' @param n number of observations
##' @param tmax maximum follow-up time
##' @param epsilon controls the censoring rate
##' 
##' @return a dataframe of five columns:
##' * \code{survtime} - survival times
##' * \code{status} - final statuses
##' * \code{age} - center and scaled age
##' * \code{Sex} - gender indicator
##' * \code{Black} - race indicator
##' @importFrom stats quantile rexp rbinom rnorm runif
##' @importFrom utils tail 
##' @export

genDataWeibull <- function(weibA, weibB, beta, n, tmax, epsilon) {
  beta <- do.call(matrix, args = list(beta, length(beta)))
  rtmax <- runif(n, 0, tmax)
  vepsilon <- rbinom(n, 1, epsilon)
  rtmax <- tmax * vepsilon + (1 - vepsilon) * rtmax
  p <- length(beta)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  x3 <- rbinom(n, 1, 0.1)
  x <- cbind(x1, x2, x3)
  U <- runif(n, 0, 1)
  tfail <- ((-log(U) / weibB * exp(-x %*% beta)) ^ (1 / weibA))
  status <- (tfail <= rtmax + 0)
  time <- tfail * status + rtmax * (1 - status)
  return(data.frame(
    survtime = time,
    status = status,
    age = x1,
    Sex = x2,
    Black = x3
  ))
}
