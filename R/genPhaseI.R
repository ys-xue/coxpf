##' Generate Phase I Profiles and Their Summaries
##' 
##' @usage genPhaseI(nprofile, baseline, alphalevel)
##' 
##' @param nprofile Number of Phase I profiles
##' @param baseline \code{"exp"} for exponential baseline hazard, \code{"weib"} for Weibull
##' baseline hazard
##' @param alphalevel overall alpha level for Phase I in determining the thresholds,
##' corresponding to \eqn{\alpha_0} in Equation 13 in the manuscript.
##' 
##' @return a list of items:
##' \itemize{
##' \item alpha - the respective alpha levels chosen for each of the four
##' statistics 
##' \item d.upper - the control limit for vertical shift
##' \item T1.upper - the control limit for max local shape deviation
##' \item T2.upper - the control limit for cumulative shape deviation
##' \item t2.upper - the control limit for Hotelling's \eqn{T^2}; if `nprofile < 100`,
##' it is the empirical quantile for Phase I \eqn{T^2} statistics; otherwise it
##' is the quantile of \eqn{\chi^2_3} distribution corresponding to alpha level
##' \item betaK - the weighted average beta vector from Phase I
##' \item muhatdelta - the median of Phase I profile centers
##' \item shatdelta - the median absolute deviation of Phase I profile centers
##' \item refShape - the estimated reference shape curve
##' \item refVariability - the estimated reference variability curve
##' \item refCor - the reference correlation matrix of Phase I covariates
##' }
##' @importFrom purrr map map_dfc map2 map2_dfc
##' @importFrom magrittr %>%
##' @importFrom stats quantile qchisq
##' @importFrom utils tail
##' @export

genPhaseI <- function(nprofile, baseline, alphalevel){
  newdata <- data.frame(survtime = 60, status = 0, age = 0, Sex = 0, Black = 0)
  
  if (baseline == "exp") {
    blockList <- map(1:nprofile, ~ genData(0.02, c(0.7, -0.5, 0.4),
                                           sample(2000:4000, 1), 60,
                                           runif(1, 0.1, 0.9)))
  } else if (baseline == "weib") {
    blockList <- map(1:nprofile, ~ genDataWeibull(0.7, 0.06, c(0.7, -0.5, 0.4), 
                                                  sample(2000:4000, 1),
                                                  60, runif(1, 0.1, 0.9)
                                                  ))
  }
  
  blockFits <- map(blockList,
                          ~coxph(Surv(survtime, status) ~ age + Sex + Black,
                                 data = .x))
  profileList <- map(blockFits, ~getProfile(.x, 1:60, newdata))
  blockbetas <- map_dfc(profileList, "beta.orig")
  blockbasehazs <- t(map_dfc(profileList, "bh.orig"))
  blockbhsds <- t(map_dfc(profileList, "bh.std"))
  betavarList <- map(profileList, "beta.var")
  betainvvarList <- betavarList %>% map(solve)
  vhList <- map(profileList, "beta.var") %>% map(~sqrt(diag(.x)))
  RList <- map(betavarList, .f = function(.x){
    diag(1/sqrt(diag(.x))) %*% .x %*% diag(1/sqrt(diag(.x)))
  })
  
  sum1 <- Reduce("+", betainvvarList)
  sum2 <- Reduce("+", map2(blockbetas, betainvvarList, `%*%`))
  betaK <- solve(sum1, t(sum2))
  
  avgCorMat <- Reduce("+", RList) / nprofile
  betadiff <- sweep(blockbetas, 1, betaK)
  betadiff <- map2_dfc(betadiff, vhList, `/`)
  t2 <- unname(diag(t(betadiff) %*% solve(avgCorMat, as.matrix(betadiff))))
  
  # ## construct a matrix of normative profiles
  deltas <- matrixStats::rowMedians(as.matrix(blockbasehazs))
  el <- sweep(blockbasehazs, 1, deltas) / blockbhsds
  ## LAD estimate of location, shape and variability
  ## search for the best value among 500 values between the max and min at each
  ## time point 
  if (baseline == "exp") {
    laddd <- ladStats(el, 0.3, 0.8, 500, 500, 1:60)
  } else if (baseline == "weib") {
    laddd <- ladStats(el, 0.7, 0.3, 500, 500, 1:60)
  }
  T1 <- unname(laddd$T1)
  T2 <- unname(laddd$T2)
  muhatdelta <- median(deltas)
  shatdelta <- mad(deltas)
  d <- abs(deltas - muhatdelta) / shatdelta
  
  ## determine $\alpha$ using the four statistics
  ## type = 8 give median unbiased quantiles

  good.alpha <- c()
  for (alpha in seq(0.001, 0.1, by = 0.00005)) {
    numvio <- length(unique(
      c(which(d > quantile(d, 1 - alpha, na.rm = TRUE, type = 8)), 
        which(T1 > quantile(T1, 1 - alpha, na.rm = TRUE, type = 8)),
        which(T2 > quantile(T2, 1 - alpha, na.rm = TRUE, type = 8)),
               which(t2 > qchisq(1 - alpha, 3)))))
    if (numvio <= ceiling(nprofile * alphalevel)) {
      good.alpha <- c(good.alpha, alpha)
    }
  }
  
  
  ## control limits:
  
  alpha <- tail(good.alpha, 1)
  
  d.upper  <- unname(quantile(d, 1 - alpha, na.rm = TRUE, type = 8))
  T1.upper <- unname(quantile(T1, 1 - alpha, na.rm = TRUE, type = 8))
  T2.upper <- unname(quantile(T2, 1 - alpha, na.rm = TRUE, type = 8))
  t2.upper <- unname(qchisq(1 - alpha, 3))
  
  return(list(alpha = alpha,
              d.upper = d.upper,
              T1.upper = T1.upper,
              T2.upper = T2.upper,
              t2.upper = t2.upper,
              betaK = betaK,
              muhatdelta = muhatdelta,
              shatdelta = shatdelta,
              refShape = laddd$shape,
              refVariability = laddd$var,
              refCor = avgCorMat))
}

