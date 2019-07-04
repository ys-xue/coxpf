# coxpf

## Overview

`coxpf` contains the code supplement for the manuscript 
**Simultaneous Monitoring for Regression Coeﬃcients and Baseline Hazard Proﬁle in
Cox Modeling of Time-to-Event Data**. 

## Installtion

```r
# use devtools to install:
# install.packages("devtools")
devtools::install_github("ys-xue/coxpf")
```

## Example code

For demonstration purpose, we present the code for replicating results in Table 1,
under the null hypothesis when `K=100` for exponential baseline, and when `K=200`
for Weibull baseline.


The output is a 5x500 matrix, with each column representing the empirical size for
one replicate. The first row denote the overall type I error, and the second to
fifth row indicate the proportion of rejections raised by each of the four statistics.


Calling `rowMeans()` on each matrix will give the average overall size, which
correspond to 0.072 and 0.062 in the first row of Table 1, as well as the average
proportion of rejections raised by each of the four statistics, corresponding
to the first and fifth rows in Supplemental Table 1.

```r
need.packages <- function(pkgs, repos = getOption("repos"), ...)
{
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new.pkgs) > 0) {
    if (is.null(repos) || repos == "@CRAN@") {
      repos <- "https://cloud.r-project.org"
    }
    install.packages(pkgs = new.pkgs, repos = repos, ...)
  }
  sapply(pkgs, function(a) {
    suppressMessages(require(a, character.only = TRUE))
  })
  invisible()
}

need.packages(c("future", "survival", "tidyverse", "matrixStats", "Rcpp",
                "RcppArmadillo", "furrr", "CoxPF"))
```

One replicate of simulation is wrapped in a `sim()` function.

```r
sim <- function(nprofile, bhtype, seed){
  set.seed(seed)
  temp <- do.call(genPhaseI, list(nprofile, bhtype, 0.05))
  d.upper <- temp[[2]]
  T1.upper <- temp[[3]]
  T2.upper <- temp[[4]]
  t2.upper <- temp[[5]]
  betaK <- temp[[6]]
  muhatdelta <- temp[[7]]
  shatdelta <- temp[[8]]
  shapes <- temp[[9]]
  variabilities <- temp[[10]]
  avgCor <- temp[[11]]
  
  if (bhtype == "exp") {
    blockList <- map(1:500, ~ genData(0.02, c(0.7, -0.5, 0.4),
                                      sample(2000:4000, 1), 60,
                                      runif(1, 0.1, 0.9)))
  } else if (bhtype == "weib") {
    blockList <- map(1:500, ~ genDataWeibull(0.7, 0.06, c(0.7, -0.5, 0.4),
                                             sample(2000:4000, 1), 60,
                                             runif(1, 0.1, 0.9)))
  }
  newdata <- data.frame(survtime = 60, status = 0, age = 0, Sex = 0, Black = 0)
  blockFits <- map(blockList,
                   ~coxph(Surv(survtime, status) ~ age + Sex + Black, data =
                            .x))
  profileList <- map(blockFits, ~getProfile(.x, 1:60, newdata))
  phase2bhs <- t(map_dfc(profileList, "bh.orig"))
  phase2betas <- map_dfc(profileList, "beta.orig")
  phase2bhsds <- t(map_dfc(profileList, "bh.std"))

  vhlist <- map(profileList, "beta.var") %>% map(~sqrt(diag(.x)))
  
  betadiff2 <- sweep(phase2betas, 1, betaK)
  betadiff2 <- map2_dfc(betadiff2, vhlist, `/`)
  HT1 <- unname(diag(t(betadiff2) %*% solve(avgCor, as.matrix(betadiff2))))
  p2medians <- matrixStats::rowMedians(phase2bhs)
  D1 <- abs(p2medians - muhatdelta) / shatdelta
  Phase2centered <- sweep(phase2bhs, 1, p2medians)
  Phase2centered <- Phase2centered / phase2bhsds
  # subtract shape
  Phase2centered <- sweep(Phase2centered, 2, shapes)
  # divide by variability
  Phase2centered <-  t(apply(Phase2centered, 1, function(x) x / variabilities))
  
  T11 <- apply(abs(Phase2centered), 1, max)
  T21 <- apply(abs(Phase2centered), 1, sum)
  violate.ind <- (D1 > d.upper |  T11 > T1.upper |  T21 > T2.upper |
                    HT1 > t2.upper)
  
  return(c(rejrate = mean(violate.ind), d1flag = mean(D1 > d.upper),
           T11flag = mean(T11 > T1.upper), T21flag = mean(T21 > T2.upper),
           HT1flag = mean(HT1 > t2.upper)))
}
```

We run the simulation studies for 500 replicates. Due to minor differences
in `R` version, the results may be slightly different. 

```r
nsims <- 500

## use parallel processing for time sake
library(furrr)
plan(multiprocess(workers = parallel::detectCores() - 1))

## replicate simulation results for null, exponential, with K=100
result100 <- future_map_dfc(1:nsims, ~sim(100, "exp", .x + 100),
                            .progress = TRUE)

round(rowMeans(result100), 3)
## results obtained using R 3.5.1 on Red Hat Linux
## [1] 0.072 0.021 0.021 0.021 0.018
## Results might be slightly different by R version
## 
## For example, sample() is a base function that can perform differently
## by R version:
## 
## On R 3.6.0, Docker image with Linux kernel
## > set.seed(1)
## > sample(2000:4000, 1)
## [1] 3016
## 
## On R 3.5.1, Docker image with Linux kernel
## > set.seed(1)
## > sample(2000:4000,1)
## [1] 2531


## replicate simulation results for null, Weibull, with K=200
result200 <- future_map_dfc(1:nsims, ~sim(200, "weib", .x + 100),
                            .progress = TRUE)


round(rowMeans(result200), 3)
## [1] 0.062 0.018 0.018 0.018 0.016


future:::ClusterRegistry("stop")
```
