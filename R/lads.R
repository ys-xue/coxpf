##' Obtain the LAD Estimations of Reference Shape and Variability
##' 
##' This is a wrapper funciton to obtain the least absolute deviation estimates
##' of reference shape and variability when we have a number of measurements,
##' each of which consist of (time/location, value) pairs.
##' 
##' @usage lads(profileMat, bw1, bw2, pc1, pc2, timevec)
##' @param profileMat A matrix or data.frame, each line of which contain one
##' set of measurements on the time grid
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
##' @return A list of three objects
##' * \code{locations} - the medians for each profile
##' * \code{shapes} - the reference shape
##' * \code{var} - the reference variability
##' 
##' @import Rcpp
##' @importFrom matrixStats rowMedians
##' @export

lads <- function(profileMat, bw1, bw2, pc1, pc2, timevec){
    ## get medians
    profileMedians <- rowMedians(profileMat)
    ## sweep profile medians from each row of matrix
    profileMat <- sweep(profileMat, 1, c(profileMedians))
    ## get shapes
    profileShapes <- do.call(lad_shape, list(profileMat, bw1, pc1, timevec))
    ## sweep profile shpes from each column of matrix
    profileMat <- sweep(profileMat, 2, c(profileShapes))
    ## get variability
    profileVars <- do.call(lad_variability, list(profileMat, bw2, pc2, timevec))
    return(list(locations = c(profileMedians), shapes = c(profileShapes),
                var = c(profileVars)))
}
