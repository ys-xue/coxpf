##' Compute Three Nonparametric Descriptive Statistics
##' 
##' Given a matrix of profiles, this function fits the LAD analysis, computes
##' the reference shape and variation curves, and returns the profile centers
##' and the matrix of relative shift. It works on top of \code{\link{lads}}.
##' 
##' @usage ladStats(profileMat, bw1, bw2, pc1, pc2, timevec)
##' 
##' @param profileMat A matrix or data.frame, each line of which contain one
##' set of measurements on the time grid, after adjusting for pointwise
##' variation
##' @param bw1 Bandwidth for RBF kernel used in estimating the reference
##' shape
##' @param bw2 Bandwidth for RBF kernel used in estimating the reference
##' variability
##' @param pc1 Precision of the grid to use in estimating the reference
##' shape
##' @param pc2 Precision of the grid to use in estimating the reference
##' variability
##' @param timevec The vector of time points where the reference shape and
##' variability curves are to be estimated
##' 
##' @return A list:
##' \itemize{
##' \item \code{D}  - profile centers
##' \item \code{T1} - local shape deviation
##' \item \code{T2} - cumulative shape deviation
##' \item \code{shape} - the reference shape
##' \item \code{variability} - the reference variability
##' }
##' @seealso \code{\link{lads}}
##' @importFrom stats mad median
##' @export
##' 


ladStats <- function(profileMat, bw1, bw2, pc1, pc2, timevec){
  estim <- lads(profileMat, bw1, bw2, pc1, pc2, timevec)
  locations <- estim$locations
  D1 <- abs(locations - median(locations)) / mad(locations)
  el1 <- sweep(profileMat, 1, locations)
  el1temp <- sweep(el1, 2, estim$shapes)
  el2 <- t(apply(el1temp, 1, function(x) x / estim$var))
  T1 <- unlist(apply(el2, 1, function(x) max(abs(x))))
  T2 <- unlist(apply(el2, 1, function(x) sum(abs(x))))
  
  return(list(D = D1, T1 = T1, T2 = T2, shape = estim$shapes,
              variability = estim$var))
}