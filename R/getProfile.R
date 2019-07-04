##' Obtain the Cumulative Baseline Hazard for A Given Time Grid
##' 
##' This function takes a \code{coxph} object, and extracts the cumulative baseline
##' hazard.
##' 
##' @usage getProfile(coxFit, times, newdata)
##' @param coxFit A \code{coxph} object
##' @param times The vector of time points where the cumulative baseline hazards
##' are to be estimated
##' @param newdata A new dataset containing the observation whose cumulative
##' baseline hazard we want to predict. If this argument is not given, an
##' hypothetical observation whose covariates equal the average of covariates
##' for all observations in the original data will be used.
##' 
##' @return A list of four objects:
##' * \code{bh.orig} -  the estimated cumulative baseline hazard
##' * \code{bh.std} - the pointwise standard deviation of \code{bh.orig}
##' * \code{beta.orig} - the coefficient vector from the Cox regression
##' * \code{beta.var} - the variance-covariance matrix of estimated coefficients
##' @import survival
##' @importFrom stats approx
##' @export


getProfile <- function(coxFit, times, newdata) {
    if (is.null(newdata)){
        sFit <- do.call(survfit, list(coxFit))
    } else {
        sFit <- survfit(coxFit, newdata = newdata)
    }
    cumhaz <- do.call(approx, list(sFit$time, -log(sFit$surv), times))
    stdhaz <- do.call(approx, list(sFit$time, sFit$std.err, times))
    
    return(list(
        bh.orig = cumhaz$y,
        bh.std = stdhaz$y,
        beta.orig = coxFit$coefficients,
        beta.var = coxFit$var
    ))
}
